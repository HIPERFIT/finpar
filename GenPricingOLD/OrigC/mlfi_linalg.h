/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mlfi_linalg.h 41326 2012-01-03 17:30:52Z jmeber $ */

#ifndef _MLFI_LINALG_H_
#define _MLFI_LINALG_H_
#include "mt19937ar.h"

#include "mlfi_pricing.h"

typedef double* vector;
typedef double* matrix;

enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};

int cholesky(matrix a, const int n);

double mlfi_ugaussian_Pinv(const double p);
void mlfi_ugaussian_Pinv_vector(vector p, const int n);

typedef struct
{
  unsigned int l;
  unsigned int *li;
  unsigned int *ri;
  unsigned int *bi;
  double *lw;
  double *rw;
  double *sd;
} mlfi_brownianbridge;

mlfi_brownianbridge *mlfi_brownianbridge_create(const double *dates, const int len);

int mlfi_brownianbridge_wiener_path(const mlfi_brownianbridge *b, enum CBLAS_TRANSPOSE trans, const matrix z, const matrix w,
				    const int rz, const int cz, const int rw, const int cw);

void mlfi_brownianbridge_normal_deviates(enum CBLAS_TRANSPOSE trans, double *w, const int r, const int c);

void mlfi_brownianbridge_free(mlfi_brownianbridge *b);

void mlfi_print_vector(const vector v, const int r);
void mlfi_print_matrix(const matrix m, const int r, const int c);

typedef struct spline {
  int n;
  int last_i;
  vector x, y, yxx, fc; /* a curve represented as a spline */
} spline;


double mlfi_spline_f (const spline *s, const double xi);
double mlfi_spline_df(const spline *s, const double xi);
double mlfi_spline_fi(const spline *s, const double xi);

void tridiag_lu_init(int n, double *al, double *ml);
void tridiag_lu_solve(int n, double *al, double *w, double *v, double *q);
void tridiag_lu_init_solve (const int n, double *a, const double *v, double *x);

void interior_scheme_1d (const int n, double *ml, const double *a_xx, const double *a_x, const double *a, const double dx, const double dt);


void vectorize_change_of_variable_1d (const double ft,
                                      const double dft,
                                      const double kt,
                                      const double dkt,
                                      const double* gx,
                                      const double* dgx,
                                      const double* d2gx,
                                      double* s,
                                      double* j_x,
                                      double* j_t,
                                      double* h_x,
                                      const int n);
void apply_change_of_variable_1d (double* a_xx, double* a_x, const double* j_x, const double* j_t, const double* h_x , const int n);
void apply_change_of_variable_2d (double* b_xx,
                                  double* b_xy,
                                  double* b_yy,
                                  double* b_x,
                                  double* b_y,
                                  const double* j_x,
                                  const double* jx_t,
                                  const double* h_x,
                                  const double* j_y,
                                  const double* jy_t,
                                  const double* h_y,
                                  const int nx,
                                  const int ny);

void tridiag_gemv(enum CBLAS_TRANSPOSE v_trans,
                  const int asize,
                  double *al,
                  const double alpha,
                  double *xbase,
                  const double beta,
                  double *ybase);
#endif
