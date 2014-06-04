/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mlfi_pricing.h 41326 2012-01-03 17:30:52Z jmeber $ */

/* MLFi pricer code C interface. This file is supposed to contain only type definitions
   and macros, no C code.

   We include everything that is necessary for generated code, so the former needs only to include
   this file */
#ifndef _MLFI_PRICING_H
#define _MLFI_PRICING_H

#include <stddef.h>
#include <float.h>
#include <math.h>

/* Keep synchronized with Mlfi_isdatypes.currency */
typedef enum {
  USD, JPY, EUR, GBP, CHF, HKD, CAD,
  AUD, DKK, NOK, SEK, CZK, IEP, MYR,
  NZD, SGD, THB, ZAR, FIM, PTE, IDR,
  TWD, EEK, HUF, ARS, BRL, CLP, FRF,
  DEM, ITL, ESP, BEF, ATS, NLG, LUF,
  GRD, ILS, KRW, LBP, MXP, PHP, PLZ,
  RUB, SAR, SKK, TRL, CNY, INR, MXN,
  TRY, PLN, IRR, AED, VEF, COP, EGP,
  NGN, PKR, RON, DZD, PEN, KZT, UAH,
  KWD, QAR, BDT, VND, MAD
} mlfi_currency;

#define MLFI_MIN_DATE 3600L
#define MLFI_MAX_DATE 168307199L

#define V_unit (mlfi_variant)0
#define V_bool(x) (mlfi_variant)(x ? 2 : 1)
#define V_int(x) (mlfi_variant)3,(mlfi_variant)x
#define V_float(x,y) (mlfi_variant)4,(mlfi_variant)x,(mlfi_variant)y
#define V_string(x) (mlfi_variant)5,(mlfi_variant)x
#define V_date(x) (mlfi_variant)6,(mlfi_variant)x
#define V_cur(x) (mlfi_variant)7,(mlfi_variant)x
#define V_tuple(x) (mlfi_variant)8,(mlfi_variant)x
#define V_list(x) (mlfi_variant)9,(mlfi_variant)x
#define V_array(x) (mlfi_variant)10,(mlfi_variant)x
#define V_None (mlfi_variant)11
#define V_Some(x) (mlfi_variant)12,(mlfi_variant)x
#define V_record(x) (mlfi_variant)13,(mlfi_variant)x
#define V_constructor(x,y) (mlfi_variant)14,(mlfi_variant)x,(mlfi_variant)y
#define V_field(x) (mlfi_variant)x

typedef enum {mlfi_false, mlfi_true} mlfi_bool;

/* lower 16 bits contain the currency, higher 16 bit the rating. */
typedef unsigned long mlfi_rated_currency;

typedef void* mlfi_variant;
typedef unsigned long mlfi_contract_index;
typedef unsigned long mlfi_slot;
typedef long mlfi_date;

/* short forms (used when generating huge files). */
#define Vr(x) V_record(x)
#define Vcc(x) V_constructor(x,0)
#define Vc(x,y) V_constructor(x,y)
#define Vf(x) V_field(x)
#define Vs(x) V_string(x)
#define mf static mlfi_variant

#ifdef _MSC_EXTENSIONS
#define MLFI_EXPORT extern __declspec(dllexport)
#define inline __inline
inline double fmax(const double x, const double y) {return x > y ? x : y;}
inline double fmin(const double x, const double y) {return x > y ? y : x;}
#else
#define MLFI_EXPORT
#define fmax(a,b) ({double _a = (a), _b = (b); _a > _b ? _a : _b; })
#define fmin(a,b) ({double _a = (a), _b = (b); _a > _b ? _b : _a; })
#endif

#define plus_day_periods(x, y) ((x) + (mlfi_date)((y) * 1440.0))
#define sub_date(x, y) ((double)((x) - (y)) / 1440.0)

/* external library routines */
double erfc(const double x);
long array_index(const int n, const int elm, const int* arr);
void order(int*, const int, const int*, const double*);

typedef struct {
  mlfi_date t;
  mlfi_date zc_t;
  mlfi_rated_currency cur;
} mlfi_cash_flow;

typedef unsigned long mlfi_pay_date_index;
typedef unsigned long mlfi_state_date_index;

/* Definitions for the new Monte-Carlo architecture. */

typedef enum {mlfi_error, mlfi_warning, mlfi_message} error_warning_message;

typedef struct {
  char *arg;
  mlfi_date t;
} mlfi_nmc_notify_parameter;

typedef struct {
  char *arg;
}  mlfi_nmc_custom_pricer;

typedef struct {
  char *arg;
  mlfi_date t;
} mlfi_nmc_deterministic_pricer;

typedef struct mlfi_nmc_model {
  void *model_data;  /* "global" data owned by the model */
  void *pricing_data; /* "global" data owned by the pricing code */
  double **underlyings;
  void (*notify_cash_flow) (const struct mlfi_nmc_model *, const mlfi_contract_index, const double, const mlfi_pay_date_index);
  void (*notify) (const struct mlfi_nmc_model *, const double, const long);
  void (*ensure_date) (const struct mlfi_nmc_model *, const mlfi_state_date_index);
  double (*custom_pricer) (const struct mlfi_nmc_model *, const mlfi_slot, const mlfi_state_date_index);
  double* deterministic_values;
  void (*error_warning_message) (const struct mlfi_nmc_model *, const error_warning_message, const char *);
} mlfi_nmc_model;

/* self-description of the pricing code, used to initialize the model */
typedef struct {
  long nb_contracts;

  long nb_state_dates;
  mlfi_date* state_dates;

  long nb_cash_flows;
  mlfi_cash_flow* cash_flows;

  long nb_std_underlyings;
  char** std_underlyings;

  long nb_notify_parameters;
  mlfi_nmc_notify_parameter* notify_parameters;

  long nb_custom_pricers;
  mlfi_nmc_custom_pricer* custom_pricers;

  long nb_deterministic_pricers;
  mlfi_nmc_deterministic_pricer* deterministic_pricers;

  long only_initial_date_pricing;

  long size_of_pricing_data;
  /* the model must allocate that many bytes for its 'pricing_data' field */

  mlfi_date max_valuation_date;
  /* the pricing code does not want to run after that date */

  void (*init) (const mlfi_nmc_model *);
  void (*trajectory) (const mlfi_nmc_model *);
}  mlfi_nmc_pricing_code;

/* Definitions for new PDE generator */

#define D_UNDERLYING 0
#define D_PATHDEP 1
#define D_BPATHDEP 2

typedef long mlfi_dimension;

typedef struct pde_dimension_info {
  long dimension_type;
  char * dimension_name;
} pde_dimension_info;

typedef struct pde_dimension_list {
  long nb_dimension;
  long * dimensions;
} pde_dimension_list;

typedef struct pde_custom_pricer_arg {
  char * arg;
} pde_custom_pricer_arg;

typedef struct pde_model {
  void* model_data; // Set by make_pde_model (custom implementation)
  mlfi_date valuation_date;
  double* (*vector_of_dimension)                (const struct pde_model*, const mlfi_date, const long);
  long    (*index_of_value)                     (const struct pde_model*, const mlfi_date, const long, const double);
  long    (*dimension_size)                     (const struct pde_model*, const long);
  void    (*discount)                           (const struct pde_model*, const long, const long, double*, const mlfi_date, const mlfi_date);
  void    (*notify_result)                      (const struct pde_model*, const long, const long, double*);
  double* (*malloc_vector)                      (const struct pde_model*, const long, const long);
  void    (*free_vector)                        (const struct pde_model*, const long, const long, double*);
  double* (*custom_pricer)                      (const struct pde_model*, const mlfi_date, const long);
  void    (*interpolation_of_value)             (const struct pde_model*, const mlfi_date, const long, const double, long*, long*, double*, double*);

 // Path-dependency support
  double  (*get_lower_bound)                    (const struct pde_model*, const mlfi_date, const long); // Lower bound of an underlying dimension on a given date
  double  (*get_upper_bound)                    (const struct pde_model*, const mlfi_date, const long); // Upper bound of an underlying dimension on a given date
  double  (*get_custom_pricer_lower_bound)      (const struct pde_model*, const mlfi_date, const long); // Lower bound of an custom pricer on a given date
  double  (*get_custom_pricer_upper_bound)      (const struct pde_model*, const mlfi_date, const long); // Upper bound of an custom pricer on a given date
  long    (*dimension_size_pathdep)             (const struct pde_model*, const long); // Size of a path-dependent dimension (typically a function that only depends on the model parameters)

  // Set by make_pde_model (default implementation)
  void*   (*malloc)                             (const struct pde_model*, const long);
  void    (*free)                               (const struct pde_model*, void*);
  void    (*discount_n)                         (const struct pde_model*, const long, double**, const long, const mlfi_date, const mlfi_date);
  void    (*free_n)                             (const struct pde_model*, long id, double**, const long);
  double**(*malloc_vector_n)                    (const struct pde_model*, const long, const long);
  double* (*vector_of_dimension_pathdep)        (const struct pde_model*, const double*, const double*, const long);
  long    (*index_of_value_pathdep)             (const struct pde_model*, const double*, const double*, const long, const double);
  void    (*interpolation_of_value_pathdep)     (const struct pde_model*, const double*, const double*, const long, const double, long *, long *, double *, double *);
  double  (*full_interpolation_of_value_pathdep)(const struct pde_model*, const double*, const double*, const long, const double, double **, long);
  void    (*printf)                             (const char *format, ...);
  void    (*flush)                              ();

} pde_model;

struct pde_model_info {
  long nb_dimensions;
  pde_dimension_info * dimensions;

  long nb_vector_infos;
  pde_dimension_list * vector_infos;
  long latest_discount_date;
  long nb_custom_pricer_args;
  pde_custom_pricer_arg * custom_pricer_args;
  long nb_discounts;
  long nb_contracts;
  void (*pde_pricer)(struct pde_model *);
};

typedef void (* pde_model_pricer) (const struct pde_model*);

typedef struct pde_model_info mlfi_pde_pricing_code;

#endif /* _MLFI_PRICING_H */
