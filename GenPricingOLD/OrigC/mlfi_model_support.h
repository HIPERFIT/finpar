/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mlfi_model_support.h 41326 2012-01-03 17:30:52Z jmeber $ */

/* This file implements types and functions that are shared amongst
   several C implementations of models. */

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

int write(int fd, const void *buf, int count);

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <time.h>
#include <sys/timeb.h>

#include "mlfi_pricing.h"
#include "mlfi_variant.h"
#include "mlfi_linalg.h"
#include "mlfi_sobol.h"
#include "mt19937ar.h"

double mlfi_erfc(double);

mlfi_currency mlfi_currency_of_string(char *s);
variant *currency(mlfi_currency cur);

/* Memory allocation */

#define mlfi_malloc(size, ty) ((size) == 0 ? NULL : (ty*) malloc((size) * sizeof(ty)))

double *mlfi_malloc_vector(size_t size, double init);
double **mlfi_malloc_matrix(size_t size1, size_t size2, double init);


/* Random number generators, Normal distribution */

typedef enum {Mersenne, Sobol} rng_t;

rng_t rng_of_string(char *s);

#define _rng(v) rng_of_string(_constructor(v)->name)

#define string_of_rng(rng) (rng == Mersenne ? "Mersenne" : "Sobol")

typedef struct {
  rng_t kind;
  union {
    qrng *sobol;
    rng mersenne;
  };
  int nb_dates;
  int dim;
} mlfi_rng;


mlfi_rng *mlfi_create_rng(rng_t k, int nb_dates, int dim);

/* fill a matrix with random normal deviates */
void mlfi_genmatrix_normal(mlfi_rng *rng, matrix zd);

/* Spline */

spline *read_spline(variant *v);

typedef struct {
  int n;
  double *x_values;
  double *y_values;
} spread;

spread* read_spread(variant*, mlfi_date);
double spread_f(spread*, double);

#define SPLINE_CURVE 0
#define AFFINE_FORWARD_CURVE 1

typedef struct {
  int curve_kind;

  /* SPLINE_CURVE case */
  spline* spline_curve;

  /* AFFINE_CURVE case */
  int nb_points;
  double *x_values;
  double *y_values;
} curve;

curve* read_curve(variant_constructor*, mlfi_date);
double curve_f(const curve*, double);

/* Step functions */
typedef struct {
  double start_date;
  double end_date;
  double value;
} step;

typedef struct {
  int length;
  step *steps;
} step_function;

step_function *read_step_function(variant*, mlfi_date);
double value_step_function(step_function* function, double t);
double integrate_step_function(step_function* function, double from, double to);
double integrate_square_step_function(step_function* function, double from, double to);

/* Driver */

#define EXsimple_price(cur,p,s) constructor("EXsimple_price", currency(cur), float(p), string(s))
#define EXvariant(v,s) constructor("EXvariant", mk_variant(v), string(s))
#define EXerror(s) constructor("EXerror", string(s))

void write_variant_to_file(char *filename, variant *result);

void write_result_exit(variant *result);

void write_error_exit(char *msg);

variant *read_command_file(char *filename);

mlfi_nmc_pricing_code *load_pricing_code(char *filename);

void run_driver(variant *(*f)(mlfi_nmc_pricing_code *, variant *), int argc, char *argv[]);


/* Callback structure */

void std_message(const struct mlfi_nmc_model * model, const error_warning_message ewn, const char *msg);

double std_custom_pricer(const mlfi_nmc_model *model, const mlfi_slot slot, const mlfi_state_date_index index);

void std_ensure_date(const mlfi_nmc_model* model, const mlfi_state_date_index index);

mlfi_nmc_model *make_model(void *md, double **u, void (*notify_cash_flow) (const struct mlfi_nmc_model *, const mlfi_contract_index, const double, const mlfi_pay_date_index), const mlfi_nmc_pricing_code *pc);


/* Monte-Carlo loop for a simple model */

void feedback(const char* s, const long length);

void run_pricing_index(mlfi_nmc_model **models, long nmodels,
                       mlfi_nmc_pricing_code *pc,
                       void (*new_trajectory)(void *),
                       void *data,
                       long niters, char *name, long * index);

void run_pricing(mlfi_nmc_model **models, long nmodels,
                 mlfi_nmc_pricing_code *pc,
                 void (*new_trajectory)(void *),
                 void *data,
                 long niters, char *name);

/* Helper for model initialization */

vector compute_deterministic_values(double (*f)(void*,variant*,mlfi_date), void* data, const mlfi_nmc_pricing_code *pc, mlfi_date today);

#define years_since(t,today) (assert(today <= t), (t - today) / (365.0 * 1440.0))
#define adapt_date(t,today) ((t) == MLFI_MIN_DATE ? today : (t))

mlfi_date *adapt_state_dates(mlfi_date today, mlfi_date *dates, long n);

mlfi_currency get_single_currency(const mlfi_nmc_pricing_code *pc);

/* For model with all payments in the same currency and a deterministic spot rate curve given as a spline,
   compute the deterministic discounts for all the payment dates */
vector discounts_single_currency(const mlfi_nmc_pricing_code *pc, mlfi_date today, curve *yc, spread* spread);

/*  Timing */

#ifdef _MSC_EXTENSIONS
typedef struct __timeb32 mlfi_timeb;
#define mlfi_ftime _ftime32
#else
typedef struct timeb mlfi_timeb;
#define mlfi_ftime ftime
#endif

#define mlfi_diff_time(t1,t2) \
  (t1.time - t2.time) * 1000 + (t1.millitm - t2.millitm)


/* Interleave sets of dates */

typedef struct {
  long ndates;
  mlfi_date *dates;
} datetbl;

datetbl *date_array(mlfi_date *dates, long k);
datetbl *single_date(mlfi_date date);

datetbl *merge_dates(datetbl *t1, datetbl *t2);

long *date_mapping(datetbl *t1, datetbl *t2);

datetbl *get_state_dates(const mlfi_date today, const mlfi_nmc_pricing_code *pc);

datetbl *get_pay_dates(mlfi_date today, mlfi_nmc_pricing_code *pc);

datetbl *_datetbl(variant *v);

datetbl *empty_datetbl();

double *year_fractions(datetbl *t, mlfi_date today);

/* Driver for pde models */

pde_model *make_pde_model(void *md);

mlfi_pde_pricing_code *load_pde_pricing_code(char *filename);

void run_pde_driver(variant *(*f)(mlfi_pde_pricing_code *, variant *), int argc, char *argv[]);

#ifdef __cplusplus
}
#endif
