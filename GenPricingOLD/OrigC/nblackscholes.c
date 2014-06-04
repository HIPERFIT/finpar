/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: nblackscholes.c 41326 2012-01-03 17:30:52Z jmeber $ */

/* A BlackScholes MonteCarlo model, that includes the following features:

   - multi-dimensional;
   - discrete dividends;
   - term-structured volatility (described by a step function);
   - interest rate described by a constant rate, a spline, or an affine forward rate function;
   - quanto adjustment;
   - non-stochastic closed formula for zero-coupons and forward rates.

   In addition to single or multiple pricing (one or several
   contracts, one or several sets of parameters), this implementation
   also supports VaR analysis (a first diffusion is applied to generate
   new sets of parameters at a later date; volatility shocks are also available).
*/

#include "mlfi_model_support.h"
#include <string.h>

#include <stdio.h>

#define get(v,field,ty) ty(_field(v,field))

/* Pre-computed data for piece-wise linear payoff custom pricers */

typedef struct{
  double a;
  double b;
  double alpha;
  double beta;
} plp_segments;

typedef struct{
  int dim;
  int pos_underlying;
  mlfi_date maturity;
  double t;
  double r_a_t;
  double r_m_t;
  double vol_t;
  plp_segments *seg;
} plp_data;

double phi(double x)
{
  return 1 - 0.5 * mlfi_erfc(x/(sqrt(2)));
}

struct driver;
typedef struct {
  /* Data read or computed by load_instance: */
  struct driver *md;
  matrix relative_dividends;
  vector spots;
  vector quanto_adjustments;
  matrix vols;
  matrix c;   /* cholesky factorization of correlation mtx */
  curve *forward_yc;  /* the spot rate curve */
  spread *spread;
  curve **asset_yc; /* spot rate curve for each asset */
  step_function **dividend_yc;
  plp_data **custom_plp;
  /* Data computed by initialize_pricing: */
  matrix drift;

  matrix trajectory;   /* 1 sample path that includes state dates +
                          dividend dates. Line (-1) contains the spots. */
  double **underlyings; /* Only state dates. */

  vector discounts;    /* multipliers for present (today) value   */
  vector vhat;         /* n * the risk-neutral expected DCF    */
  step_function **volfunction;
  mlfi_nmc_model *model; /* the callback structure bound to this instance */
  /* Used during instance initialization: */
} instance;

typedef struct driver {
  FILE* debug_file;
  mlfi_nmc_pricing_code* mi;
  int mean;
  unsigned long nb_paths;
  mlfi_date today;
  int dim;
  int nb_path_dates;
  mlfi_currency cur;

  /* Initialized by initialize_model */
  mlfi_rng *rng;
  mlfi_brownianbridge *bb;

  /* Data local to each trajectory, but shared amongts instances.
     Pre-allocated in initialize_model. */
  matrix zd;           /* normal deviates before brownian bridge    */
  matrix z;            /* normal deviates after brownian bridge     */

  int nb_instances;
  instance **instances;

  mlfi_nmc_model **models; /* the callback structures of all instances */

  /* Used during instance initialization: */
  int nb_dividend_dates;
  int nb_state_dates;
  long *dividend_date_map; /* map dividend dates to path dates */
  long *state_date_map; /* map state dates to path dates */
  vector fdates;    /* dates as floats in years with today=0.0 */
  double start_date;
} driver;


double deterministic_pricers(void *data, variant* v, mlfi_date t) {
  variant_constructor *kind = _constructor(v);
  instance *inst = data;
  if (strcmp(kind->name,"Zcb_risky") == 0) {
    mlfi_date mat;
    double amount, taut, taumat, x, xmat;

    mat = _date(kind->arguments[1]);
    amount = _float(kind->arguments[2]);

    taut   = years_since(t, inst->md->today);
    taumat = years_since(mat, inst->md->today);

    x    = taut   * (curve_f(inst->forward_yc, taut) + spread_f(inst->spread, taut));
    xmat = taumat * (curve_f(inst->forward_yc, taumat) + spread_f(inst->spread, taumat));

    return amount * exp(- (xmat - x));
  } else if (strcmp(kind->name,"Zcb_risk_free") == 0) {
    mlfi_date mat;
    double amount, taut, taumat, x, xmat;

    mat = _date(kind->arguments[1]);
    amount = _float(kind->arguments[2]);

    taut   = years_since(t, inst->md->today);
    taumat = years_since(mat, inst->md->today);

    x    = taut   * (curve_f(inst->forward_yc, taut));
    xmat = taumat * (curve_f(inst->forward_yc, taumat));

    return amount * exp(- (xmat - x));
  } else if (strcmp(kind->name,"Forward") == 0) {
    int t1, t2;
    double daycfraction, xt1, xt2, taut1, taut2;

    t1 = _date(kind->arguments[1]);
    t2 = _date(kind->arguments[2]);
    daycfraction = _float(kind->arguments[3]);

    taut1  = years_since(t1, inst->md->today);
    taut2  = years_since(t2, inst->md->today);

    xt1  = taut1 * curve_f(inst->forward_yc, taut1);
    xt2  = taut2 * curve_f(inst->forward_yc, taut2);

    return (exp(xt2 - xt1) - 1.0) / daycfraction;
  } else {
    assert(0 && "Unknown deterministic pricer");
    return 0.0;
  }
}

void initialize_custom_plp(const instance *inst, const mlfi_slot slot, const mlfi_nmc_pricing_code *pc) {
  int i;
  double t;
  variant *arg;
  variant_constructor *kind;
  mlfi_date today = inst->md->today;
  arg = mlfi_parse_variant(pc->custom_pricers[slot].arg);
  assert(arg && "Cannot parse custom pricer string");
  kind = _constructor(arg);
  if (strcmp(kind->name,"PLP") == 0) {
    plp_data *pdata = mlfi_malloc(1,plp_data);
    variant *params = kind->arguments[0];
    variant *segments = get(params, segments, _list);
    int pos_s;
    char *name = get(params,underlying,_string);
    for (pos_s = 0; pos_s < pc->nb_std_underlyings; pos_s++)
      if (strcmp(pc->std_underlyings[pos_s], name) == 0)
        break;
    assert((pos_s < pc->nb_std_underlyings) && "Underlying not found");
    pdata->pos_underlying = pos_s;
    pdata->maturity = get(params, maturity, _date);
    t = years_since(pdata->maturity, today);
    pdata->t = t;
    pdata->seg = mlfi_malloc(segments->len,plp_segments);
    pdata->dim = segments->len;
    pdata->r_a_t = t * (curve_f(inst->asset_yc[pos_s],t) + value_step_function(inst->dividend_yc[pos_s],t));
    pdata->r_m_t = t * (curve_f(inst->forward_yc, t) + spread_f(inst->spread, t));
    pdata->vol_t = integrate_square_step_function(inst->volfunction[pos_s], today, pdata->maturity);
    for(i=0; i < segments->len; i++) {
      variant *vi = _nth(segments, i);
      pdata->seg[i].a = get(vi, a, _float);
      pdata->seg[i].b = get(vi, b, _float);
      pdata->seg[i].alpha = get(vi, alpha, _float);
      pdata->seg[i].beta = get(vi, beta, _float);
    }
    inst->custom_plp[slot]= pdata;
  }
  else assert("error not an available custom pricer");
}

void notify_cash_flow(const mlfi_nmc_model *model,
                      const mlfi_contract_index contract_number,
                      const double amount,
                      const mlfi_pay_date_index date_index)
{
  instance *inst = (instance*) model->model_data;
  FILE* debug_file = inst->md->debug_file;
  inst->vhat[contract_number] += amount * inst->discounts[date_index];
  if (NULL != debug_file) {
    mlfi_nmc_pricing_code* mi = inst->md->mi;
    mlfi_date today = inst->md->today;
    double date = years_since(adapt_date(mi->cash_flows[date_index].t, today), today);
    fprintf(debug_file, "CashFlow, %3li, %12g, %12g, %12g\n", contract_number, amount, date, inst->discounts[date_index]);
  }
}

void new_trajectory(void *data)
{
  driver *md = data;
  int dim = md->dim;
  int npathdates = md->nb_path_dates;
  int i, j, k, l, m;
  double *z = md->z;
  FILE* debug_file = md->debug_file;

  if (dim==0 || npathdates == 0) return;
  /* so we can function as a 'deterministic' model */

  /* This is shared by all instances */
  mlfi_genmatrix_normal(md->rng, md->zd);

  mlfi_brownianbridge_wiener_path(md->bb, CblasNoTrans, md->zd, z, npathdates, dim, npathdates, dim);

  for (i = dim * npathdates - 1; i >= dim; i--) z[i] -= z[i-dim];

  /* Instance-specific processing */
  for (m = 0; m < md->nb_instances; m++) {
    instance *inst = md->instances[m];
    for (i = 0; i < npathdates; i++)  {
      for (j = 0; j < dim; j++) {
        double temp = 0.;
        k = dim * i + j;
        for (l = 0; l <= j; l++) temp += inst->c[dim * j + l] * z[dim * i + l];
        temp = exp(temp * inst->vols[k] + inst->drift[k]);
        inst->trajectory[k] = inst->trajectory[k - dim] * temp;
      }
    }
  }

  if (NULL != debug_file) {
    fprintf(debug_file, "\nNewTrajectory\n\n");
    for (m = 0; m < md->nb_instances; m++) {
      instance *inst = md->instances[m];
      for (i = 0; i < npathdates; i++)  {
        for (j = 0; j < dim; j++) {
          k = dim * i + j;
          fprintf(debug_file, "Underlying, %25s, %12g, %12g\n",
                  md->mi->std_underlyings[j],
                  md->fdates[i],
                  inst->trajectory[k]);
        }
      }
    }
    fprintf(debug_file, "\n");
  }
}

/* Load the data from the command file and initialize model */

instance *load_instance(const variant *v, driver *md, const mlfi_nmc_pricing_code *mi) {
  int i,j;
  variant *assets;
  matrix correlation;
  int dim = md->dim;
  int npathdates = md->nb_path_dates;
  instance *inst = mlfi_malloc(1, instance);
  inst->md = md;
  /* Read the spot rate curve */
  inst->forward_yc = read_curve(get(v,curve,_constructor), md->today);
  inst->spread = read_spread(get(v,spread,_list), md->today);
  assets = get(v,assets,_array); assert(assets->len == dim);
  if (!dim) return inst;
  /* Asset data */
  inst->spots = mlfi_malloc(dim, double);
  inst->quanto_adjustments = mlfi_malloc(dim, double);
  inst->vols = mlfi_malloc(dim * npathdates, double);
  inst->relative_dividends = mlfi_malloc_vector(dim * npathdates, 1.);
  inst->asset_yc = mlfi_malloc(dim,curve*);
  inst->dividend_yc = mlfi_malloc(dim,step_function*);
  inst->volfunction = mlfi_malloc(dim,step_function*);

  /* Correlation */
  correlation = mlfi_malloc(dim * dim, double);

  for (i=0; i < dim; i++) {
    step_function *vol;
    variant *dividends;
    variant *asset = assets->vlist[i];
    int correl_len;
    double* correl;

    assert(strcmp(get(asset,name,_string), mi->std_underlyings[i]) == 0);
    vol = read_step_function(get(asset, vol, _record), md->today);
    inst->volfunction[i] = vol;
    inst->spots[i] = get(asset, spot, _float);

    correlation[i*dim+i] = 1.0;
    correl_len = i;
    correl = _float_array(_field(asset, correl), &correl_len);
    for (j=0; j < i; j++){
      correlation[i*dim+j] = correlation[j*dim+i] = correl[j];
    }

    inst->quanto_adjustments[i] = get(asset, quanto_adjustment, _float);
    for (j=0; j < npathdates; j++){
      if (0 == j) {
        double var = 0;
        if (0 != md->fdates[j])
          var = integrate_square_step_function(vol, 0, md->fdates[j]) / md->fdates[j];
        inst->vols[dim * j + i] = sqrt(var);
      }
      else {
        double var = integrate_square_step_function(vol, md->fdates[j-1], md->fdates[j]) / (md->fdates[j] - md->fdates[j-1]);
        inst->vols[dim * j + i] = sqrt(var);
      }
     }

    dividends = get(asset,dividends,_array);
    assert(dividends->len == md->nb_dividend_dates);
    for (j=0; j < md->nb_dividend_dates; j++) {
      double div = _float(dividends->varray[j]);
      long k = md->dividend_date_map[j]*dim+i;
      inst->relative_dividends[k] = 1.0 - div;
    }

    inst->asset_yc[i] = read_curve(get(asset,asset_curve,_constructor), md->today);
    inst->dividend_yc[i] = read_step_function(get(asset,continuous_dividend_curve,_record), md->today);
  }

  if (cholesky(correlation, dim) != 0) {write_error_exit("Cholesky decomposition failed");};
  inst->c = correlation;

  inst->custom_plp = mlfi_malloc(mi->nb_custom_pricers, plp_data*);
  for(i=0; i < mi->nb_custom_pricers; i++)
    initialize_custom_plp(inst, i, mi);

  return inst;
}

void compute_blackscholes_drift(instance *inst) {
  driver *md = inst->md;
  int npathdates = md->nb_path_dates;
  int dim = md->dim;
  double t0 = md->start_date;
  vector bm;
  int i,j;

  inst->drift = mlfi_malloc(npathdates * dim, double);
  bm = mlfi_malloc_vector(dim, 0.);
  if (t0 != 0)
    for (j=0; j < dim; j++) {
      double vol;
      vol = integrate_step_function(inst->volfunction[j], 0, t0);
      bm[j] = t0 * (curve_f(inst->asset_yc[j], t0) + value_step_function(inst->dividend_yc[j], t0)) - inst->quanto_adjustments[j] * vol;
    }
  for (i=0; i < npathdates; i++) {
    double t = md->fdates[i] + md->start_date;
    double dt = i==0 ? md->fdates[i] : (md->fdates[i]-md->fdates[i-1]);
    for (j=0; j < dim; j++) {
      double vol = integrate_step_function(inst->volfunction[j], 0, t);
      double f = t * (curve_f(inst->asset_yc[j], t) + value_step_function(inst->dividend_yc[j], t)) - inst->quanto_adjustments[j] * vol;
      int k = dim * i + j;
      inst->drift[k] =
        f - bm[j]
        - dt *
        (0.5 * inst->vols[k] * inst->vols[k])
        + log(inst->relative_dividends[k]);
      bm[j] = f;
    }
  }
  free(bm);
}

void initialize_diffusion_instance(instance *inst) {
  int dim = inst->md->dim;
  compute_blackscholes_drift(inst);
  inst->trajectory = mlfi_malloc((1 + inst->md->nb_path_dates) * dim, double);
  /* The first 'line' (will be -1) is a copy of spots. This is to avoid
     a test in the inner loop of new_trajectory */
  memcpy(inst->trajectory, inst->spots, sizeof(double) * dim);
  inst->trajectory += dim;
}

void initialize_diffusion_driver(driver *md, const rng_t rng, const int extra_rn) {
  /* Pre-allocate memory for new_trajectory */
  md->z = mlfi_malloc(md->nb_path_dates * md->dim, double);
  md->zd = mlfi_malloc((md->nb_path_dates + extra_rn) * md->dim, double);

  /* Create random number generator */
  md->rng = mlfi_create_rng(rng, md->nb_path_dates + extra_rn, md->dim);
  md->bb = mlfi_brownianbridge_create(md->fdates, md->nb_path_dates);
}

double custom_pricer(const mlfi_nmc_model *model, const mlfi_slot slot, const mlfi_state_date_index index) {
  double t,r_m,s,r_a,vol2,variance;
  double res = 0;
  int pos_s;
  instance *inst = (instance*)model->model_data;
  double t0 = inst->md->fdates[index];
  int i;
  plp_data *temp;
  temp = inst->custom_plp[slot];
  t = temp->t - t0;
  r_m = (temp->r_m_t - t0 * (curve_f(inst->forward_yc, t0) + spread_f(inst->spread, t0)))/t;
  pos_s = temp->pos_underlying;
  r_a = ((temp->r_a_t - t0 * (curve_f(inst->asset_yc[pos_s], t0) + value_step_function(inst->dividend_yc[pos_s], t0)))/t) ;
  s = inst->underlyings[index][pos_s];
  vol2 = integrate_square_step_function(inst->volfunction[pos_s], t0, temp->t);
  variance = sqrt(vol2);
  for(i=0; i < temp->dim; i++) {
    double z0,z1,u0,u1;
    double a=temp->seg[i].a;
    double b=temp->seg[i].b;
    double alpha=temp->seg[i].alpha;
    double beta=temp->seg[i].beta;
    z1 = (log((s*exp(r_a*t))/b) + 0.5*vol2)/variance;
    u1 = (log((s*exp(r_a*t))/b) - 0.5*vol2)/variance;
    if (a != 0){
      z0 = (log((s*exp(r_a*t))/a) + 0.5*vol2)/variance;
      u0 = (log((s*exp(r_a*t))/a) - 0.5*vol2)/variance;
      res += exp(-1*r_m*t)* ((exp(r_a *t)*alpha*s*(phi(z0)  - phi(z1)))+(beta*(phi(u0)  - phi(u1))));
    }
    else{
      res += exp(-1*r_m*t)* ((exp(r_a *t)*alpha*s*(1 - phi(z1)))+(beta*(1 - phi(u1))));
    }
  }
  return res;
}



void initialize_pricing(instance *inst, const mlfi_nmc_pricing_code *mi, const int start_date, const double *past) {
  int i,j;
  driver *md = inst->md;
  int dim = md->dim;
  double *past_vals;
  /* Discounts and discounted payments */
  inst->vhat = mlfi_malloc_vector(mi->nb_contracts, 0.);
  inst->discounts = discounts_single_currency(mi, md->today, inst->forward_yc, inst->spread);

  inst->underlyings = NULL; /* needed for codegen in the deterministic case (cf codegen_trajectory_method) */
  if (!dim) goto deterministic_model;
  /* allow us to function as a 'deterministic' model */

  initialize_diffusion_instance(inst);

  j = 0;
  for (i = 0; i < md->nb_state_dates; i++)
    if (md->state_date_map[i] < start_date) j++; else break;
  past_vals = mlfi_malloc(dim * j, double);

  inst->underlyings = mlfi_malloc(md->nb_state_dates, double*);
  for (i = 0; i < md->nb_state_dates; i++) {
    int j = md->state_date_map[i];
    if (j < start_date) {
      memcpy(&past_vals[i*dim], &past[j*dim], sizeof(double) * dim);
      inst->underlyings[i] = &past_vals[i*dim];
    }
    else {
    	assert( (i==(j-start_date)) && "Trajectory does NOT match underlyings!!!\n\n"); // Cosmin
    	inst->underlyings[i] = &inst->trajectory[(j - start_date)*dim];
    }
  }

 deterministic_model:

  /* Create the callback structure and run deterministic pricers */
  inst->model = make_model(inst, inst->underlyings, notify_cash_flow, mi);
  inst->model->deterministic_values =
    compute_deterministic_values(&deterministic_pricers,inst,mi,md->today);
  inst->model->custom_pricer = &custom_pricer;
}



driver *initialize_model(const variant *v, mlfi_nmc_pricing_code *mi)
{

  driver *md = mlfi_malloc(1, driver);
  int i, j, k;
  int dim = mi->nb_std_underlyings;
  int npathdates;
  double *fdates;
  datetbl *dividend_dates, *state_dates, *path_dates, *var_dates;
  long *var_dates_map;
  mlfi_date today;
  variant_constructor *kind;
  rng_t rng;
  char* debug_file = get(v, debug_file, _string);
  if (*debug_file != 0)
    md->debug_file = fopen(debug_file, "w");
  else
    md->debug_file = NULL;
  md->mi = mi;
  md->dim = dim;
  md->start_date = 0;
  md->cur = get_single_currency(mi);
  md->nb_paths = get(v,nb_paths,_int);
  today = md->today = get(v,today,_date);
  assert("max_valuation_date" && md->today <= mi->max_valuation_date);
  kind = get(v,instances,_constructor);

  rng = get(v,random_number_generator,_rng);

  dividend_dates = get(v,dividend_dates,_datetbl);
  md->nb_dividend_dates = dividend_dates->ndates;

  state_dates = get_state_dates(today, mi);
  md->nb_state_dates = state_dates->ndates;

  if (strcmp(kind->name,"VaR_analysis") == 0) {
    mlfi_date t = get(kind->arguments[0],var_date,_date);
    if (t < today)
      write_error_exit("VaR analysis at a date before pricing date");
    var_dates = single_date(t);
  }
  else
    var_dates = empty_datetbl();

  /* Compute path dates (interleave state dates and dividend dates) */
  path_dates = merge_dates(merge_dates(state_dates, dividend_dates), var_dates);
  md->dividend_date_map = date_mapping(dividend_dates, path_dates);
  md->state_date_map = date_mapping(state_dates, path_dates);
  var_dates_map = date_mapping(var_dates, path_dates);
  md->nb_path_dates = npathdates = path_dates->ndates;
  md->fdates = fdates = year_fractions(path_dates, today);

  /* Read all instances */
  if (strcmp(kind->name,"Explicit") == 0) {
    variant *instances = kind->arguments[0];
    md->nb_instances = instances->len;
    md->instances = mlfi_malloc(instances->len, instance*);

    for (i=0; i < instances->len; i++) {
      md->instances[i] = load_instance(_nth(instances,i),md,mi);
      initialize_pricing(md->instances[i],mi,0,NULL);
    }
    md->mean = 0;
  } else if (strcmp(kind->name,"VaR_analysis") == 0) {
    /* Use the current driver "md" so that load_instance will
       use the real number of path dates and the real dividend_date_map. */

    instance *instsim =
      load_instance(get(kind->arguments[0],var_init,_record),md,mi);

    int nb_simul = get(kind->arguments[0],var_nb_simul,_int);
    long nbpd = var_dates_map[0] + 1;
    driver *mdsim = mlfi_malloc(1, driver);
    variant *vol_shock = get(kind->arguments[0],var_vol_shock,_list);
    int shock = 0;

    if (vol_shock->len != 0) {
      assert(vol_shock->len == dim);
      for (j = 0; j < dim; j++)
        if (_float(vol_shock->vlist[j]) != 0) shock = 1;
    }

    mdsim->today = md->today;
    mdsim->debug_file = NULL;
    mdsim->dim = md->dim;
    mdsim->fdates = fdates;
    mdsim->nb_instances = 1;
    mdsim->instances = &instsim;
    mdsim->nb_path_dates = nbpd;
    mdsim->start_date = 0;
    instsim->md = mdsim;
    initialize_diffusion_driver(mdsim, rng, 1);
    initialize_diffusion_instance(instsim);

    md->nb_instances = nb_simul;
    md->instances = mlfi_malloc(nb_simul, instance*);
    md->nb_path_dates = npathdates - nbpd;
    md->fdates = fdates + nbpd;
    md->start_date = fdates[nbpd-1];
    for (i = 0; i < md->nb_path_dates; i++)  md->fdates[i] -= md->start_date;

    for (i = 0; i < nb_simul; i++) {
      instance *inst = mlfi_malloc(1, instance);
      new_trajectory((void*) mdsim);
      inst->md = md;
      inst->relative_dividends = instsim->relative_dividends + nbpd * dim;
      inst->spots = &instsim->trajectory[(nbpd - 1) * dim];
      inst->volfunction = instsim->volfunction;
      inst->quanto_adjustments = instsim->quanto_adjustments;
      if (!shock)
        inst->vols = instsim->vols + nbpd * dim;
      else {
        inst->vols = mlfi_malloc(md->nb_path_dates * dim, double);
        for (j = 0; j < dim; j++) {
          double shock = _float(vol_shock->vlist[j]) * mdsim->zd[mdsim->nb_path_dates * dim + j];
          for (k = 0; k < md->nb_path_dates; k++)
            inst->vols[k * dim + j] = instsim->vols[(nbpd + k) * dim + j] + shock;
          /* TODO: protect against vol < 0 ? */
        }
      }
      inst->c = instsim->c;
      inst->spread = instsim->spread;
      inst->forward_yc = instsim->forward_yc;
      inst->asset_yc = instsim->asset_yc;
      inst->dividend_yc = instsim->dividend_yc;
      initialize_pricing(inst,mi,nbpd,instsim->trajectory);
      md->instances[i] = inst;
    }
    md->mean = 0;
    //md->mean = 1;
  } else {
    write_error_exit("Pricing kind not supported");
  }

  initialize_diffusion_driver(md,rng, 0);

  md->models = mlfi_malloc(md->nb_instances, mlfi_nmc_model*);
  for (i=0; i < md->nb_instances; i++) md->models[i] = md->instances[i]->model;
  return md;
}

variant *blackscholes(mlfi_nmc_pricing_code *pc, variant *cmd) {
  driver *md;
  variant *result;
  int i,j,k,niters;
  mlfi_timeb t1,t2;
  long int elapsed;
  assert(pc->nb_state_dates > 0);
  mlfi_ftime(&t1);
  md = initialize_model(cmd, pc);
  assert(pc->init);
  assert(pc->trajectory);
  /* run the simulation */
  niters = (pc->only_initial_date_pricing ? 1 : md->nb_paths);
  run_pricing(md->models, md->nb_instances, pc,
              new_trajectory, md,
              niters, "monte-carlo pricing");
  /* produce the result */
  mlfi_ftime(&t2);
  elapsed = mlfi_diff_time(t2,t1);
  if (NULL != md->debug_file) fclose(md->debug_file);

  if (!md->mean) {
    result = create_list(pc->nb_contracts * md->nb_instances);
    k = 0;
    for (i=0; i < pc->nb_contracts; i++) {
      for (j = 0; j < md->nb_instances; j++) {
        result->vlist[k++] =
          EXsimple_price(md->cur,md->instances[j]->vhat[i]/niters,
                         bprintf("MC BS (C) %ddim, %s, %d paths, %ld ms",
                                 pc->nb_std_underlyings,
                                 string_of_rng(md->rng->kind),
                                 md->nb_paths, elapsed));
      }
    }
  } else {
    double sum = 0;
    result = create_list(1);
    for (j = 0; j < md->nb_instances; j++)  sum += md->instances[j]->vhat[0];
    result->vlist[0] =
      EXsimple_price(md->cur,sum / md->nb_instances / niters,
                     bprintf("VAR MC BS (C) %ddim, %s, %d paths, %ld ms",
                             pc->nb_std_underlyings,
                             string_of_rng(md->rng->kind),
                             md->nb_paths, elapsed));
  }
  return result;
}

int main(int argc, char *argv[]) {
  run_driver(blackscholes,argc,argv);
  return 1; // Cosmin
}
