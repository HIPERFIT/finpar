/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mlfi_model_support.c 42117 2012-01-26 16:26:34Z mmouraret $ */

/* This file implements types and functions that are shared amongst
   several C implementations of models. */

#include "mlfi_model_support.h"
#include "mlfi_erfc.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <time.h>
#include <sys/timeb.h>

/* Keep synchronized with Mlfi_isdatypes.currency */
static char *mlfi_currencies[] = {
  "USD", "JPY", "EUR", "GBP", "CHF", "KHD", "CAD",
  "AUD", "DKK", "NOK", "SEK", "CZK", "IEP", "MYR",
  "NZD", "SGD", "THB", "ZAR", "FIM", "PTE", "IDR",
  "TWD", "EEK", "HUF", "ARS", "BRL", "CLP", "FRF",
  "DEM", "ITL", "ESP", "BEF", "ATS", "NLG", "LUF",
  "GRD", "ILS", "KRW", "LBP", "MXP", "PHP", "PLZ"
  "RUB", "SAR", "SKK", "TRL", "CNY", "INR", "MXN",
  "TRY", "PLN", "IRR", "AED", "VEF", "COP", "EGP",
  "NGN", "PKR", "RON", "DZD", "PEN", "KZT", "UAH",
  "KWD", "QAR", "BDT", "VND", "MAD", "XAU", "XAG",
  "XPT", "XPD",
  0
};

#include "mlfi_variant.c"

mlfi_currency mlfi_currency_of_string(char *s) {
  int i;
  for (i = 0; mlfi_currencies[i]; i++)
    if (strncmp(s,mlfi_currencies[i],3) == 0) return i;
  return (-1);
}

variant* currency(mlfi_currency cur) {
  return constructor(mlfi_currencies[cur]);
}


/* Memory allocation */

double *mlfi_malloc_vector(size_t size, double init) {
  double *v = mlfi_malloc(size, double);
  int i;
  for (i = 0; i < size; i++) v[i] = init;
  return v;
}

double **mlfi_malloc_matrix(size_t size1, size_t size2, double init) {
  double **v = mlfi_malloc(size1, double*);
  int i;
  for (i = 0; i < size1; i++) v[i] = mlfi_malloc_vector(size2, init);
  return v;
}


/* Random number generators */

rng_t rng_of_string(char *s) {
  if (strcmp(s,"Mersenne") == 0) return Mersenne;
  else if (strcmp(s,"Sobol") == 0) return Sobol;
  else { printf("rng must be == Mersenne or Sobol, not %s\n", s); exit(2); }
}

mlfi_rng *mlfi_create_rng(rng_t k, int nb_dates, int dim) {
  mlfi_rng *rng = mlfi_malloc(1, mlfi_rng);
  rng->kind = k;
  rng->nb_dates = nb_dates;
  rng->dim = dim;

  if (k == Sobol) {
    rng->sobol = qrng_alloc(qrng_sobol(), dim*nb_dates);
  }
  else {
    unsigned long init_key[] = {217182818L}; /* note this seed matches the default seed in MLFi */
    mt_init_by_array(&rng->mersenne, &init_key[0], 1);
  }
  return rng;
}

/* fill a matrix with random (0,1) uniform deviates */
void mlfi_genmatrix_uniform(mlfi_rng *rng, matrix zd) {
  if (rng->kind == Sobol)
    qrng_get(rng->sobol, zd);
  else {
    long i,j;
    long dim = rng->dim;
    for (i=0; i < rng->nb_dates; i++)
      for (j=0; j < dim; j++)
	zd[i*dim+j] = mt_genrand_real3(&rng->mersenne);
  }
}

/* fill a matrix with random normal deviates */
void mlfi_genmatrix_normal(mlfi_rng *rng, matrix zd) {
  mlfi_genmatrix_uniform(rng, zd);
  mlfi_ugaussian_Pinv_vector(zd, rng->nb_dates * rng->dim);
}

/* Spline */

spline *read_spline(variant *v) {
  spline *s = mlfi_malloc(1, spline);
  s->last_i = 0;
  s->n = (-1);
  s->x = _float_array(_field(v,x), &s->n);
  s->y = _float_array(_field(v,y), &s->n);
  s->yxx = _float_array(_field(v,y''), &s->n);
  s->fc = _float_array(_field(v,fc), &s->n);
  return s;
}

spread *read_spread(variant *list, mlfi_date today) {
  spread *result = mlfi_malloc(1, spread);
  int i;
  long length = list->len;
  result->n = length;
  result->x_values = mlfi_malloc(length, double);
  result->y_values = mlfi_malloc(length, double);
  for (i=0; i<length; i++) {
    variant *current_step = _nth(list, i);
    result->x_values[i] = years_since(_date(_fst(current_step)), today);
    result->y_values[i] = _float(_snd(current_step));
  }
  return result;
}

double spread_f(spread *s, double t) {
  long low = -1;
  long up = s->n;
  double *x = s->x_values;
  double *y = s->y_values;
  while (low+1 < up) {
    long mid = (low + up) / 2;
    if (x[mid] <= t) low = mid;
    else up = mid;
  }
  if (0 == up) {
    return y[0];
  }
  else if (s->n == up) {
    return y[s->n-1];
  }
  else {
    double result = y[low] * (x[up]-t)/(x[up]-x[low]) + y[up] * (t-x[low])/(x[up]-x[low]);
    return result;
  }
}

curve* read_curve(variant_constructor *kind, mlfi_date today) {
  curve *c = mlfi_malloc(1, curve);
  if (0==strcmp(kind->name, "Spline")) {
    c->curve_kind = SPLINE_CURVE;
    c->spline_curve = read_spline(kind->arguments[0]);
  }
  else if (0==strcmp(kind->name, "Affine")) {
    variant *list = kind->arguments[0];
    int i;
    long length = list->len;
    c->curve_kind = AFFINE_FORWARD_CURVE;
    c->nb_points = length;
    c->x_values = mlfi_malloc(length, double);
    c->y_values = mlfi_malloc(length, double);
    for (i=0; i<length; i++) {
      variant *current_step = _nth(list, i);
      c->x_values[i] = years_since(_date(_fst(current_step)), today);
      c->y_values[i] = _float(_snd(current_step));
    }
  }
  else
    assert (0 && "Unknown curve kind.");
  return c;
}

double curve_f(const curve* c, double t) {
  if (SPLINE_CURVE==c->curve_kind)
    return mlfi_spline_f(c->spline_curve, t);
  else if (AFFINE_FORWARD_CURVE==c->curve_kind) {
    long n;
    double res = 0;
    double *x_values = c->x_values;
    double *y_values = c->y_values;
    double dt = t - x_values[0];
    if (dt == 0) return y_values[0];
    for (n=0; n<c->nb_points-1; n++) {
      double height = (y_values[n]+y_values[n+1])/2;
      if (x_values[n] <= t && t < x_values[n+1]) {
        res += (t-x_values[n]) * height;
        break;
      }
      else
        res += (x_values[n+1]-x_values[n]) * height;
    }
    if (n==c->nb_points-1)
      res += (t-x_values[n])*y_values[n];
    return res/dt;
  }
  else
    assert (0 && "Unknown curve kind.");
  return 0; /* avoid warning */
}

step_function *read_step_function(variant *steps, mlfi_date today) {
  double first_step = _float(_field(steps, first_step));
  variant *v = _field(steps, next_steps);
  long i;
  variant *current_step;
  step_function *sf = mlfi_malloc(1, step_function);
  long length = (v-> len) + 1;
  sf->length = length;
  sf->steps = mlfi_malloc(length, step);
  sf->steps[0].value = first_step;
  sf->steps[0].start_date = 0.0;
  for (i = 1; i < length; i++) {
    current_step = _nth(v, i - 1);
    sf->steps[i].value = _float(_snd(current_step));
    sf->steps[i].start_date = years_since(_date(_fst(current_step)), today);
    sf->steps[i-1].end_date = sf->steps[i].start_date;
  }
  sf->steps[length-1].end_date = 300.;
  return sf;
}

double value_step_function(step_function* sf, double t) {
  long i;
  for (i = 0; i < sf->length; i++) {
    double end_date = sf->steps[i].end_date;
    double v = sf->steps[i].value;
    if (t < end_date) return v;
  }
  return sf->steps[sf->length - 1].value;
}

double integrate_step_function(step_function* sf, double from, double to) {
  long i;
  double result = 0;
  for (i = 0; i < sf->length; i++) {
    double start_date = sf->steps[i].start_date;
    double end_date = sf->steps[i].end_date;
    double v = sf->steps[i].value;
    double dcf = 0;
    if ((from < start_date) && (start_date <= to) && (to <= end_date))
      dcf = to - start_date;
    else if ((from < start_date) && (end_date < to))
      dcf = end_date - start_date;
    else if ((start_date <= from) && (to <= end_date))
      dcf = to - from;
    else if ((start_date <= from) && (from <= end_date) && (end_date < to))
      dcf = end_date - from;
    result += dcf * v;
  }
  return result;
}

double integrate_square_step_function(step_function* sf, double from, double to) {
  /* todo: factorize with the previous function. */
  long i;
  double result = 0;
  for (i = 0; i < sf->length; i++) {
    double start_date = sf->steps[i].start_date;
    double end_date = sf->steps[i].end_date;
    double v = sf->steps[i].value;
    double dcf = 0;
    if ((from < start_date) && (start_date <= to) && (to <= end_date))
      dcf = to - start_date;
    else if ((from < start_date) && (end_date < to))
      dcf = end_date - start_date;
    else if ((start_date <= from) && (to <= end_date))
      dcf = to - from;
    else if ((start_date <= from) && (from <= end_date) && (end_date < to))
      dcf = end_date - from;
    result += dcf * v * v;
  }
  return result;
}

/* Driver */

void write_variant_to_file(char *filename, variant *result) {
  FILE *file = fopen(filename, "w");
  mlfi_strbuf *buf = mlfi_strbuf_create();
  assert(file);
  mlfi_print_variant(buf, result);
  fprintf(file, "%s\n", buf->data);
  fclose(file);
}

char *result_filename = NULL;

void write_result_exit(variant *result) {
  assert(result_filename);
  write_variant_to_file(result_filename, result);
  exit(0);
}

void write_error_exit(char *msg) {
  variant *result = create_list(1);
  result->vlist[0] = EXerror(msg);
  write_result_exit(result);
}

variant *read_command_file(char *filename) {
  variant *v =
    mlfi_parse_variant
    (mlfi_strbuf_getfree(mlfi_strbuf_readfile(mlfi_strbuf_create(), filename)));
  assert(v && "Cannot parse variant in command file");
  return v;
}

mlfi_nmc_pricing_code *load_pricing_code_dll(char *filename) {
  mlfi_nmc_pricing_code* (*get_pricing_code)();

#ifdef _WIN32
  HMODULE d;
  d = LoadLibrary(filename);
  assert("LoadLibrary" && d);
  *(void **)(&get_pricing_code) = GetProcAddress(d, "get_pricing_code");
  assert("GetProcAddress (get_pricing_code)" && get_pricing_code);
#else
  void *d;
#ifdef HIPERFIT
  char  filename_dll[80];
  strcpy(filename_dll, "./");
  strcat(filename_dll, filename);
  strcat(filename_dll, ".dll");
  //char* filename_dll = strcat(filename, ".dll");
  d = dlopen(filename_dll, RTLD_NOW);
#else
  d = dlopen(filename, RTLD_NOW);
#endif
  assert("dlopen" && d);
  *(void **)(&get_pricing_code) = dlsym(d, "get_pricing_code");
  assert("dlsym get_pricing_code" && get_pricing_code);
#endif
  return get_pricing_code();
}

static char* find_sep(char *s) {
  while (*s != 0) {
    if (s[0] == '\n' && s[1] == '=' && s[2] == '=' && s[3] == '=' && s[4] == '\n') {
      s[0] = 0;
      return (s + 5);
    }
    s++;
  }
  assert("cannot find separator" && 0);
  return NULL;
}

/* Note: a lot of code duplication with mlfi_native_stubs.c related to codegen. */

#define CODEGEN_ALIGN 16
#define intnat ptrdiff_t
static char* codegen_magic32 = "LEXIFI_MONTE_CARLO_PRICING_x86_2";
static char* codegen_magic64 = "LEXIFI_MONTE_CARLO_PRICING_x64_2";
static int codegen_size = 0;
static char* codegen_alloc = NULL;
static char* codegen_code = NULL;
static char* codegen_data_alloc = NULL;
static char* codegen_data = NULL;
static mlfi_nmc_pricing_code *codegen_infos = NULL;
static double** codegen_matrix = NULL;
typedef void (*codegen_execute_fun)(double **);
typedef codegen_execute_fun (*codegen_relocate_fun)(void **,char *);
static codegen_execute_fun codegen_fun = NULL;
static mlfi_nmc_model *codegen_last_model = NULL;

static int of_hexdigit(char c) {
  if ('0' <= c  && c <= '9') return (c - '0');
  if ('a' <= c  && c <= 'f') return (c - 'a') + 10;
  if ('A' <= c  && c <= 'F') return (c - 'A') + 10;
  assert("of_hexdigit" && 0);
  exit(2);
}

static void codegen_cash_flow(intnat cidx, intnat slot, double x) {
  codegen_last_model->notify_cash_flow(codegen_last_model, cidx, x, slot);
}

static void codegen_notify(double x, intnat slot) {
  codegen_last_model->notify(codegen_last_model, x, slot);
}

static double codegen_custom(intnat slot, intnat index) {
  return codegen_last_model->custom_pricer(codegen_last_model, slot, index);
}

static void codegen_order(intnat* dest, const long n, const intnat* arr, const double* fs)
{
  intnat i,j;
  intnat a;
  double fa;
  dest[0]=arr[0];
  for (j=1;j<n;j++){
    a=arr[j];
    fa=fs[a];
    i=j-1;
    while (i >= 0 && fs[dest[i]] > fa) {
      dest[i+1]=dest[i];
      i--;
    }
    dest[i+1]=a;
  }
}

static intnat codegen_array_index(intnat n, intnat elm, intnat* arr)
{
  intnat i;
  for (i=0; i<n; i++)
    if (arr[i] == elm) return i;
  return (-1);
}

static intnat *codegen_permutate(intnat length, intnat *idxs, double *vals, intnat *result) {
  codegen_order(result, length, idxs, vals);
  return result;
}

static double codegen_power(double x, double y) {
  return pow(x, y);
}

static double codegen_log(double x) {
  return log(x);
}

static double codegen_erfc(double x) {
  return mlfi_erfc(x);
}

static double codegen_floor(double x) {
  return floor(x);
}

static void *codegen_function_table[]= {
  codegen_cash_flow,
  codegen_custom,
  codegen_notify,
  codegen_permutate,
  codegen_array_index,
  codegen_power,
  codegen_log,
  codegen_erfc,
  codegen_floor
};

static void load_codegen(char *code, int size, int size_of_pricing_data) {
  if (codegen_alloc != NULL) {
#ifdef _WIN32
    VirtualFree(codegen_alloc, codegen_size + CODEGEN_ALIGN - 1, MEM_RELEASE);
#else
    free(codegen_alloc);
#endif
    codegen_alloc = NULL;
  }

  codegen_size = size;
#ifdef _WIN32
  codegen_alloc = (char*)VirtualAlloc(NULL, codegen_size + CODEGEN_ALIGN - 1, MEM_COMMIT, PAGE_EXECUTE_READWRITE);
#else
  codegen_alloc = malloc(codegen_size + CODEGEN_ALIGN - 1);
#endif
  assert(codegen_alloc);
  codegen_code = codegen_alloc + CODEGEN_ALIGN - 1;
  codegen_code = codegen_code - ((long)codegen_code % CODEGEN_ALIGN);
  memcpy(codegen_code, code, codegen_size);


  if (codegen_data_alloc != NULL) {
    free(codegen_data_alloc);
    codegen_data_alloc = NULL;
  }

  codegen_data_alloc = malloc(size_of_pricing_data + CODEGEN_ALIGN - 1);
  assert(codegen_data_alloc);
  codegen_data = codegen_data_alloc + CODEGEN_ALIGN - 1;
  codegen_data = codegen_data - ((long)codegen_data % CODEGEN_ALIGN);

  codegen_fun = *((codegen_relocate_fun)codegen_code)(codegen_function_table, codegen_data);
}

#define get(v,field,ty) ty(_field(v,field))

static void codegen_init_method(const mlfi_nmc_model *model) {
  if (codegen_matrix) {
    free(codegen_matrix);
    codegen_matrix = NULL;
  }
  codegen_matrix = malloc(sizeof(double*) * (1 + codegen_infos->nb_state_dates));
  codegen_last_model = NULL;
}

static void codegen_trajectory_method(const mlfi_nmc_model *model) {
  int i = 0;
  if (model != codegen_last_model) {
    codegen_last_model = (mlfi_nmc_model*) model;
    codegen_matrix[0] = model->deterministic_values;
    if (model->underlyings)
      for (i = 0; i < codegen_infos->nb_state_dates; i++)
        codegen_matrix[i + 1] = model->underlyings[i];
  }
  codegen_fun(codegen_matrix);
}

mlfi_nmc_pricing_code *load_pricing_code_codegen(char *filename) {
  mlfi_nmc_pricing_code *pc;
  variant *infos;
  int codelen, i;
  char *code;
  char *data0 = mlfi_strbuf_getfree(mlfi_strbuf_readfile(mlfi_strbuf_create(), filename));
  char *data1 = find_sep(data0);
  char *data2 = find_sep(data1);
  char *data3 = find_sep(data2);
  char *magic;
  int len;

  switch (sizeof(intnat)) {
  case 4: magic = codegen_magic32; break;
  case 8: magic = codegen_magic64; break;
  default: assert("codegen not supported on this platform" && 0);
  }
  if (strcmp(data0, magic) != 0)
    assert("wrong magic number in pricing code" && 0);

  infos = mlfi_parse_variant(data1);
  assert(infos && "Cannot parse pricing infos variant in codegen file");

  pc = malloc(sizeof(mlfi_nmc_pricing_code));
  codegen_infos = pc;
  pc->nb_contracts = get(infos, nb_contracts, _int);
  pc->state_dates = _date_array_len(_field(infos, state_dates), len);
  pc->nb_state_dates = len;
  pc->size_of_pricing_data = get(infos, size_of_pricing_data, _int);
  {
    variant *arr = _field(infos, cash_flows);
    int nb = mlfi_variant_len(arr);

    pc->nb_cash_flows = nb;
    pc->cash_flows = malloc(nb * sizeof(mlfi_cash_flow));
    for (i = 0; i < nb; i++) {
      variant *elt = _nth(arr, i);
      mlfi_currency cur = mlfi_currency_of_string(_constructor(_fst(_thd(elt)))->name);
      int rating = _int(_snd(_thd(elt)));
      pc->cash_flows[i].t = _date(_fst(elt));
      pc->cash_flows[i].zc_t = _date(_snd(elt));
      pc->cash_flows[i].cur = cur + (rating << 16);
    }
  }
  pc->std_underlyings = _string_array_len(_field(infos, all_std_underlyings), len);
  pc->nb_std_underlyings = len;
  pc->only_initial_date_pricing = get(infos, only_initial_date_pricing, _bool);
  pc->max_valuation_date = get(infos, max_valuation_date, _date);

  {
    variant *arr = _field(infos, deterministic_pricers);
    int nb = mlfi_variant_len(arr);

    pc->nb_deterministic_pricers = nb;
    pc->deterministic_pricers = malloc(nb * sizeof(mlfi_nmc_deterministic_pricer));
    for (i = 0; i < nb; i++) {
      variant *elt = _nth(arr, i);
      pc->deterministic_pricers[i].t = _date(_snd(elt));
      pc->deterministic_pricers[i].arg = mlfi_string_of_variant(_variant(_fst(elt)));
    }
  }

  {
    variant *arr = _field(infos, notify_parameters);
    int nb = mlfi_variant_len(arr);

    pc->nb_notify_parameters = nb;
    pc->notify_parameters = malloc(nb * sizeof(mlfi_nmc_notify_parameter));
    for (i = 0; i < nb; i++) {
      variant *elt = _nth(arr, i);
      pc->notify_parameters[i].t = _date(_snd(elt));
      pc->notify_parameters[i].arg = mlfi_string_of_variant(_variant(_fst(elt)));
    }
  }

  {
    variant *arr = _field(infos, custom_pricers);
    int nb = mlfi_variant_len(arr);

    pc->nb_custom_pricers = nb;
    pc->custom_pricers = malloc(nb * sizeof(mlfi_nmc_custom_pricer));
    for (i = 0; i < nb; i++) {
      variant *elt = _nth(arr, i);
      pc->custom_pricers[i].arg = mlfi_string_of_variant(_variant(elt));
    }
  }

  codelen = strlen(data2) / 2;
  code = malloc(codelen);
  for (i = 0; i < codelen; i++)
    code[i] = (of_hexdigit(data2[i * 2]) << 4) + of_hexdigit(data2[i * 2 + 1]);
  load_codegen(code, codelen, pc->size_of_pricing_data);
  free(code);

  pc->init = codegen_init_method;
  pc->trajectory = codegen_trajectory_method;

  return pc;
}

mlfi_nmc_pricing_code *load_pricing_code(char *filename) {
  int len = strlen(filename);
  if (len > 8 && strcmp(filename + len - 8, ".codegen") == 0)
    return load_pricing_code_codegen(filename);
  else
    return load_pricing_code_dll(filename);
}



void run_driver(variant *(*f)(mlfi_nmc_pricing_code *, variant *), int argc, char *argv[]) {
  if (argc != 4) {
    printf ("usage: %s <result file> <library name> <command file>\n",
            argv[0]);
    exit (1);
  }
  result_filename = argv[1];
  write_result_exit(f(load_pricing_code(argv[2]), read_command_file(argv[3])));
  exit(0);
}


/* Callback structure */

void std_message(const struct mlfi_nmc_model * model, const error_warning_message ewn, const char *msg) {
  assert("error warning message" && 0);
  /* TODO */
}

double std_custom_pricer(const mlfi_nmc_model *model, const mlfi_slot slot, const mlfi_state_date_index index) {
  write_error_exit("Unknown custom pricer");
  return 0.;
}

void std_ensure_date(const mlfi_nmc_model* model, const mlfi_state_date_index index) {
  /* Do nothing */
}

void std_notify(const mlfi_nmc_model* model, const double v, const long index) {
  write_error_exit("Notification not supported");
}

mlfi_nmc_model *make_model(void *md, double **u, void (*notify_cash_flow) (const struct mlfi_nmc_model *, const mlfi_contract_index, const double, const mlfi_pay_date_index), const mlfi_nmc_pricing_code *pc) {
  mlfi_nmc_model *model = mlfi_malloc(1, mlfi_nmc_model);
  model->model_data = (void*) md;
  model->pricing_data = malloc(pc->size_of_pricing_data);
  model->underlyings = u;
  model->notify_cash_flow = notify_cash_flow;
  model->notify = std_notify;
  model->ensure_date = std_ensure_date;
  model->custom_pricer = std_custom_pricer;
  model->deterministic_values = NULL;
  model->error_warning_message = std_message;
  return model;
}

/* Monte-Carlo loop for a simple model */

void feedback(const char* s, const long length){ write(1, s, length); }


#define RUN_ORIG 0

#if RUN_ORIG == 0

/**********************************/
/**** HIPERFIT_VERSION (COSMIN) ***/
/**********************************/
#if 0
void run_pricing_indexGPU (
		mlfi_nmc_model **models, long nmodels,
        mlfi_nmc_pricing_code *pc,
        void (*new_trajectory)(void *),
        void *data,
        long niters, char *name, long * index);
#endif

void mainLoopGPU (
		mlfi_nmc_model **models, long nmodels,
		mlfi_nmc_pricing_code *pc,
		void (*new_trajectory)(void *),
		void *data,
		long niters, char *name, long * index
);


int parallel_predicate_holds(void* data);



void run_pricing_index(mlfi_nmc_model **models, long nmodels,
                 mlfi_nmc_pricing_code *pc,
                 void (*new_trajectory)(void *),
                 void *data,
                 long niters, char *name, long * index) {

  long step = (niters > 20 ? niters / 20 : 1);
  long i,j;
  *index = 0;
  feedback("\n##[20:",7);
  feedback(name,strlen(name));
  feedback("\n",1);

  //printf("Cosmin num models: %d %d\n\n\n", nmodels, niters);

  for (j = 0; j < nmodels; j++) pc->init(models[j]);

  if( parallel_predicate_holds(data) && !RUN_ORIG) { // COSMIN's branch
	  mainLoopGPU(models, nmodels, pc, new_trajectory, data, niters, name, index);
  } else {

	  /** Cosmin: this seems to be the main loop
	   *          Questions:
	   */
	  for (i = 0; i < niters; i++) {
		//if (i % step == 0) feedback("##\n",3);
		new_trajectory(data);
		for (j = 0; j < nmodels; j++) pc->trajectory(models[j]);
		(*index)++;
	  }
  }
  feedback("##]\n",4);
}

#else

void run_pricing_index(mlfi_nmc_model **models, long nmodels,
                 mlfi_nmc_pricing_code *pc,
                 void (*new_trajectory)(void *),
                 void *data,
                 long niters, char *name, long * index) {

  long step = (niters > 20 ? niters / 20 : 1);
  long i,j;
  *index = 0;
  feedback("\n##[20:",7);
  feedback(name,strlen(name));
  feedback("\n",1);
  for (j = 0; j < nmodels; j++) pc->init(models[j]);

  for (i = 0; i < niters; i++) {
    if (i % step == 0) feedback("##\n",3);
    new_trajectory(data);
    for (j = 0; j < nmodels; j++) pc->trajectory(models[j]);
    (*index)++;
  }
  feedback("##]\n",4);
}
#endif


void run_pricing(mlfi_nmc_model **models, long nmodels,
                 mlfi_nmc_pricing_code *pc,
                 void (*new_trajectory)(void *),
                 void *data,
                 long niters, char *name) {
  long index;
  run_pricing_index(models, nmodels, pc, new_trajectory, data, niters, name, &index);
}

/* Helper for model initialization */

vector compute_deterministic_values(double (*f)(void*,variant*,mlfi_date), void* data, const mlfi_nmc_pricing_code *pc, mlfi_date today) {
  vector v;
  long i;
  if (!pc->nb_deterministic_pricers) return NULL;
  v = mlfi_malloc(pc->nb_deterministic_pricers, double);
  for (i=0; i < pc->nb_deterministic_pricers; i++) {
    variant *arg = mlfi_parse_variant(pc->deterministic_pricers[i].arg);
    mlfi_date t = pc->deterministic_pricers[i].t;
    assert(v && "Cannot parse variant string");

    if (t == MLFI_MIN_DATE) t = today;
    else assert(today <= t && "deterministic pricer in the past");

    v[i] = f(data, arg, t);
  }
  return v;
}

mlfi_date *adapt_state_dates(mlfi_date today, mlfi_date *dates, long n) {
  /* Note: only the first state_date might be == MLFI_MIN_DATE */
  mlfi_date *adates = mlfi_malloc(n, mlfi_date);
  long i;
  for (i=0; i<n; i++) adates[i] = adapt_date(dates[i],today);
  return adates;
}

mlfi_currency get_single_currency(const mlfi_nmc_pricing_code *pc) {
  mlfi_currency cur;
  long i;
  assert(pc->nb_cash_flows > 0);
  cur = (mlfi_currency) pc->cash_flows[0].cur; /* CHECK: mlfi_rated_currency -> mlfi_currency? */
  for (i=1; i < pc->nb_cash_flows; i++) assert(pc->cash_flows[i].cur == cur);
  return cur;
}

/* For model with all payments in the same currency and a deterministic spot rate curve given as a spline,
   compute the deterministic discounts for all the payment dates */
vector discounts_single_currency(const mlfi_nmc_pricing_code *pc, mlfi_date today, curve *yc, spread* spread)
{
  vector discounts = mlfi_malloc(pc->nb_cash_flows, double);
  long i;
  for (i=0; i < pc->nb_cash_flows; i++) {
    double df = years_since(adapt_date(pc->cash_flows[i].t, today), today);
    double r  = curve_f(yc, df) + spread_f(spread, df);
    discounts[i] = exp(-r*df);
  }
  return discounts;
}

/* Interleave sets of dates */

datetbl *date_array(mlfi_date *dates, long k) {
  datetbl *res = mlfi_malloc(1, datetbl);
  res->dates = dates;
  res->ndates = k;
  return res;
}

datetbl *single_date(mlfi_date date) {
  mlfi_date *dates = mlfi_malloc(1, mlfi_date);
  dates[0] = date;
  return date_array(dates, 1);
}

datetbl *merge_dates(datetbl *t1, datetbl *t2) {
  long i,j,k;
  mlfi_date *dates = mlfi_malloc(t1->ndates + t2->ndates, mlfi_date);
  for (i = 0, j = 0, k = 0; i < t1->ndates || j < t2->ndates; ) {
    if ((i < t1->ndates) && (j == t2->ndates || t1->dates[i] < t2->dates[j]))
      dates[k++] = t1->dates[i++];
    else if ((j < t2->ndates) && (i == t1->ndates || t2->dates[j] < t1->dates[i]))
      dates[k++] = t2->dates[j++];
    else
      { j++; dates[k++] = t1->dates[i++]; }
  }
  return date_array(dates,k);
}

long *date_mapping(datetbl *t1, datetbl *t2) {
  long i,j;
  long *res = mlfi_malloc(t1->ndates, long);
  for (i = 0, j = 0; i < t1->ndates; i++) {
    while (t2->dates[j] < t1->dates[i]) j++;
    res[i] = j;
  }
  return res;
}

datetbl *get_state_dates(const mlfi_date today, const mlfi_nmc_pricing_code *pc) {
  return date_array(adapt_state_dates(today, pc->state_dates,
                                      pc->nb_state_dates),
                    pc->nb_state_dates);
}

datetbl *get_cash_flows(mlfi_date today, mlfi_nmc_pricing_code *pc) {
  long n = pc->nb_cash_flows;
  mlfi_date *dates = mlfi_malloc(n, mlfi_date);
  long i;
  for (i = 0; i < n; i++) dates[i] = pc->cash_flows[i].t;
  return date_array(adapt_state_dates(today, dates, n), n);
}

datetbl *_datetbl(variant *v) {
  int k;
  mlfi_date *dates = _date_array_len(v,k);
  return date_array(dates, k);
}

datetbl *empty_datetbl() {
  return date_array(NULL, 0);
}

double *year_fractions(datetbl *t, mlfi_date today) {
  long i;
  double *fdates = mlfi_malloc(t->ndates, double);
  for (i=0; i < t->ndates; i++) {
    if (t->dates[i] < today) write_error_exit("Path date before pricing date");
    fdates[i] = years_since(t->dates[i], today);
  }
  return fdates;
}

/* Driver for pde models */

mlfi_pde_pricing_code *load_pde_pricing_code(char *filename) {
  mlfi_pde_pricing_code* (*get_pricing_code)();
#ifdef _WIN32
  HMODULE d;
  d = LoadLibrary(filename);
  assert("LoadLibrary" && d);
  *(void **)(&get_pricing_code) = GetProcAddress(d, "get_pde_pricing_code");
  assert("GetProcAddress (get_pde_pricing_code)" && get_pricing_code);
#else
  void *d;
  d = dlopen(filename, RTLD_NOW);
  assert("dlopen" && d);
  *(void **)(&get_pricing_code) = dlsym(d, "get_pde_pricing_code");
  assert("dlsym get_pde_pricing_code" && get_pricing_code);
#endif
  return get_pricing_code();
}

void run_pde_driver(variant *(*f)(mlfi_pde_pricing_code *, variant *), int argc, char *argv[]) {
  if (argc != 4) {
    printf ("usage: %s <result file> <library name> <command file>\n",
            argv[0]);
    exit (1);
  }
  result_filename = argv[1];
  write_result_exit(f(load_pde_pricing_code(argv[2]), read_command_file(argv[3])));
  exit(0);
}
