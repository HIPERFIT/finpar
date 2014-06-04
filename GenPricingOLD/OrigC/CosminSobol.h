#ifndef COSMIN_SOBOL_H
#define COSMIN_SOBOL_H

#include "mlfi_model_support.h"
#include "poly.h"


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
/*
double phi(double x)
{
  return 1 - 0.5 * mlfi_erfc(x/(sqrt(2)));
}
*/
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
  curve **dividend_yc;
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


#endif
