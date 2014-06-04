/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

static char *cvsid="$Id: mlfi_linalg.c 41326 2012-01-03 17:30:52Z jmeber $";

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "mlfi_linalg.h"

/* Linear Algebra definitions which are used by both MLFi and native pricers (e.g. written in C) */

/* Cholesky */
int cholesky(matrix a, const int n){
  int i, j, k;
  double sum, *diagonal;
  diagonal = malloc(n * sizeof(double));
  for(i=0;i<n;i++) {
    for(j=i;j<n;j++) {
      for(sum=a[i*n+j],k=i-1;k>=0;k--) sum -= a[i*n+k]*a[j*n+k];
      if(i == j) {
        if(sum <= 0.0) return 1;
        diagonal[i] = sqrt(sum);
      }
      else {
        a[j*n+i] = sum/diagonal[i];
        a[i*n+j] = a[j*n+i];
      } /* end if */
    } /* end for */
  } /* end for */
  for(i=0;i<n;i++)
    a[i*n+i] = diagonal[i];
  free(diagonal);
  return 0;
}

#define rat_eval(a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7) \
    (x*(x*(x*(x*(x*(x*(x*a7+a6)+a5)+a4)+a3)+a2)+a1)+a0)/ \
    (x*(x*(x*(x*(x*(x*(x*b7+b6)+b5)+b4)+b3)+b2)+b1)+b0)

inline double small_case(const double q)
{
  double x = 0.180625 - q * q;
  return q * rat_eval(
                      3.387132872796366608,
                      133.14166789178437745,
                      1971.5909503065514427,
                      13731.693765509461125,
                      45921.953931549871457,
                      67265.770927008700853,
                      33430.575583588128105,
                      2509.0809287301226727,

                      1.0,
                      42.313330701600911252,
                      687.1870074920579083,
                      5394.1960214247511077,
                      21213.794301586595867,
                      39307.89580009271061,
                      28729.085735721942674,
                      5226.495278852854561);
}

inline double intermediate(const double r)
{
  double x = r - 1.6;
  return rat_eval(
                  1.42343711074968357734,
                  4.6303378461565452959,
                  5.7694972214606914055,
                  3.64784832476320460504,
                  1.27045825245236838258,
                  0.24178072517745061177,
                  0.0227238449892691845833,
                  7.7454501427834140764e-4,

                  1.0,
                  2.05319162663775882187,
                  1.6763848301838038494,
                  0.68976733498510000455,
                  0.14810397642748007459,
                  0.0151986665636164571966,
                  5.475938084995344946e-4,
                  1.05075007164441684324e-9);
}

inline double tail(const double r)
{
  double x = r - 5.0;
  return rat_eval(
                  6.6579046435011037772,
                  5.4637849111641143699,
                  1.7848265399172913358,
                  0.29656057182850489123,
                  0.026532189526576123093,
                  0.0012426609473880784386,
                  2.71155556874348757815e-5,
                  2.01033439929228813265e-7,

                  1.0,
                  0.59983220655588793769,
                  0.13692988092273580531,
                  0.0148753612908506148525,
                  7.868691311456132591e-4,
                  1.8463183175100546818e-5,
                  1.4215117583164458887e-7,
                  2.04426310338993978564e-15);
}

double mlfi_ugaussian_Pinv(const double p)
{
  double dp = p - 0.5;

  /* These two cases sometimes occur as limit cases in closed forms */
  if (p==1) return exp(0.0)/0.0; /* exp is MSVC workaround */
  if (p==0) return -exp(0.0)/0.0;

  if (fabs(dp) <= 0.425) return small_case(dp);

  {
    double pp = dp < 0. ? p : (1.0 - p);
    double r = sqrt (- log(pp));
    double x = r <= 5.0 ? intermediate(r) : tail(r);
    return dp < 0. ? -x : x;
  }
}

/* "vectorized" version of preceeding one, without limit case check. */
void mlfi_ugaussian_Pinv_vector(vector p, const int n)
{
  int i;
  for (i=0; i < n; i++, p++)
    {
      double dp = *p - 0.5;
      if (fabs(dp) <= 0.425)
        {
          *p = small_case(dp);
        }
      else
        {
          double pp = dp < 0. ? *p : (1.0 - *p);
          double r = sqrt (- log(pp));
          double x = r <= 5.0 ? intermediate(r) : tail(r);
          if (dp < 0.) {*p = -x;} else {*p = x;};
        }
    }
}

static int brownianbridge_create_aux(mlfi_brownianbridge *b, int n, const double *dates, const unsigned int i, const unsigned int k)
{
  unsigned int j = i + (k-i)/2;
  if (j > i && j < k) {
    double ti, tj, tk;

    b->li[n] = i;
    b->bi[n] = j;
    b->ri[n] = k;

    ti = i==0 ? 0.0 : dates[i-1];
    tj = dates[j-1];
    tk = dates[k-1];

    if (tk == ti) { /*  gracefully handle the case where adjacent dates in the list are identical */
      if (tj != tk) return -1; /* failwith("brownianbridge_create_aux: non-monotonicity in date list"); */

      b->lw[n] = 0.0;
      b->rw[n] = 0.0;
      b->sd[n] = 0.0;
    }
    else {
      b->lw[n] = (tk-tj)/(tk-ti);
      b->rw[n] = (tj-ti)/(tk-ti);
      b->sd[n] = sqrt( (tj-ti)*(tk-tj)/(tk-ti) );
    }

    n = brownianbridge_create_aux(b, n+1, dates, i, j);
    n = brownianbridge_create_aux(b, n, dates, j, k);
  }
  return n;
}

mlfi_brownianbridge *mlfi_brownianbridge_create(const double *dates, const int l)
{
  mlfi_brownianbridge *b;

  b = (mlfi_brownianbridge *) malloc(sizeof(mlfi_brownianbridge));

  if (!b) return NULL; /* failwith("brownianbridge_create: memory allocation 1"); */
  b->l = l;
  b->li = (unsigned int *) malloc(sizeof(unsigned int) * l);
  b->ri = (unsigned int *) malloc(sizeof(unsigned int) * l);
  b->bi = (unsigned int *) malloc(sizeof(unsigned int) * l);

  b->lw = (double *) malloc(sizeof(double) * l);
  b->rw = (double *) malloc(sizeof(double) * l);
  b->sd = (double *) malloc(sizeof(double) * l);

  if (!b->li || !b->ri || !b->bi || !b->lw || !b->rw || !b->sd)
    return NULL; /* failwith("brownianbridge_create: memory allocation 2"); */

  if (l > 0) {

    b->bi[0] = l;
    b->sd[0] = sqrt(dates[l-1]);
    brownianbridge_create_aux(b, 1, dates, 0, l);
  }
  return b;
}


int mlfi_brownianbridge_wiener_path(const mlfi_brownianbridge *b, enum CBLAS_TRANSPOSE trans, const matrix z, const matrix w,
				    const int rz, const int cz, const int rw, const int cw)
{

  if (trans==CblasNoTrans) {
    if (b->l != rz)
      return -1; /* number of rows doesn't match length of bridge */

    if (b->l > 0) {
      size_t m;

      for (m=0; m < cz; m++) {
	size_t i;
	w[(b->bi[0]-1) * cz + m] = b->sd[0] * z[m]; /* TODO: eliminate int mul */
	for(i=1; i < rz; i++) {
	  int j=b->li[i]-1;
	  int k=b->ri[i]-1;
	  int l=b->bi[i]-1;
	  double wk = w[k*cz+m]; /* idem */
	  double zi = z[i*cz+m]; /* idem */
	  w[l*cz+m] = /* idem */
	    j == -1
	    ?                          b->rw[i] * wk + b->sd[i] * zi
	    : b->lw[i] * w[j*cz+m] +  b->rw[i] * wk + b->sd[i] * zi; /* idem */
	}
      }
    }
  }
  else {
    if (b->l != cz)
      return -2; /* bridge size doesn't match no. columns */

    if (b->l > 0) {
      size_t m;

      for (m=0; m < rz; m++) {
	size_t i;
	w[m*cz + (b->bi[0]-1)] = b->sd[0] * z[m*cz];
	for(i=1; i < cz; i++) {
	  int j=b->li[i]-1;
	  int k=b->ri[i]-1;
	  int l=b->bi[i]-1;
	  double wk = w[m*cz+k]; /* idem */
	  double zi = z[m*cz+i]; /* idem */
	  w[m*cz+l] = /* idem */
	    j == -1
	    ?                         b->rw[i] * wk + b->sd[i] * zi
	    : b->lw[i] * w[m*cz+j] +  b->rw[i] * wk + b->sd[i] * zi; /* idem */
	}
      }
    }
  }
  return 0;
}

void mlfi_brownianbridge_normal_deviates(enum CBLAS_TRANSPOSE trans, double *w, const int r, const int c)
{
  int i, j;
  if (trans==CblasNoTrans) {
    if (r > 0)
      for (j=0; j < c; j++) {
	for (i=r-1; i > 0; i--) {
	  w[i*c+j] -= w[(i-1)*c+j]; /* idem */
	}
      }
  }
  else {
    if (c > 0)
      for (j=0; j < r; j++) {
	for (i=c-1; i > 0; i--) {
	  w[j*c+i] -= w[j*c+i-1]; /* idem */
	}
      }
  }
}

void mlfi_brownianbridge_free(mlfi_brownianbridge *b)
{
  free(b->li);
  free(b->ri);
  free(b->bi);
  free(b->lw);
  free(b->rw);
  free(b->sd);
  free(b);
}

void mlfi_print_vector(const vector v, const int r)
{
  int i;
  for (i=0; i < r; i++) {
    printf("%16.12g ", v[i]);
  }
  printf("\n");
}


void mlfi_print_matrix(const matrix m, const int r, const int c)
{
  int i, j;
  for (i=0; i < r; i++) {
    for (j=0; j < c; j++) printf("%16.12g ", m[i*c+j]);
    printf("\n");
  }
}

/* spline curve handling */

static int find(const spline *s, const double xi)
{
  int k;
  int klo=0, khi=s->n-1;

 start:
  k = (khi + klo) / 2;
  if (khi - klo > 1) {
    if (s->x[k] > xi) {
      khi=k; goto start;
    }
    else {
      klo=k; goto start;
    }
  }
  else return k;
}

double mlfi_spline_f(const spline *s, const double xi)
{
  unsigned long klo, khi;
  double h, a, b, h2, c, d;
  double r;

  klo = find(s, xi);
  khi = klo + 1L;

  h = s->x[khi] - s->x[klo];
  assert (h > 0.0);

  a = (s->x[khi] - xi) / h;
  b = (xi - s->x[klo]) / h;
  h2 = h*h;
  c = (1.0 / 6.0) * (a*a*a - a) * h2;
  d = (1.0 / 6.0) * (b*b*b - b) * h2;
  r = a*s->y[klo] + b*s->y[khi] + c*s->yxx[klo] + d*s->yxx[khi];

  return r;
}


double mlfi_spline_df(const spline *s, const double xi)
{
  unsigned long klo, khi;
  double h, a, b;
  double dadx, dbdx, dcdx, dddx;

  klo = find(s, xi);
  khi = klo + 1L;
  h = s->x[khi] - s->x[klo];
  assert (h > 0.0);

  a = (s->x[khi] - xi) / h;
  b = (xi - s->x[klo]) / h;
  dadx = -1.0/h;
  dbdx =  1.0/h;
  dcdx = h / 6.0 - h * a * a / 2;
  dddx = h * b * b / 2 - h / 6.0;

  return  dadx*s->y[klo] + dbdx*s->y[khi] + dcdx*s->yxx[klo] + dddx*s->yxx[khi];
}


static inline double ia (vector x, unsigned long klo, unsigned long khi, double tau, double h, double h2)
{  return x[khi] * tau / h - tau * tau / (2*h);}

static inline double ia3 (vector x, unsigned long klo, unsigned long khi, double tau, double h, double h2)
{  return (1.0/h) * (x[khi]*x[khi]*x[khi] * tau - (3.0/2.0)*tau*tau*x[khi]*x[khi] + tau*tau*tau*x[khi] - tau*tau*tau*tau/4.0); }


static inline double ib (vector x, unsigned long klo, unsigned long khi, double tau, double h, double h2)
{  return tau * tau / (2*h) - x[klo] * tau / h ; }


static inline double ib3 (vector x, unsigned long klo, unsigned long khi, double tau, double h, double h2)
{  return (1.0/h) * (tau*tau*tau*tau/4.0 - tau*tau*tau*x[klo] + (3.0/2.0)*tau*tau*x[klo]*x[klo] - x[klo]*x[klo]*x[klo] * tau); }


static inline double ic (vector x, unsigned long klo, unsigned long khi, double tau, double h, double h2)
{  return (1.0/6.0)*h2*(ia3(x, klo, khi, tau, h, h2) - ia(x, klo, khi, tau, h, h2));}


static inline double id (vector x, unsigned long klo, unsigned long khi, double tau, double h, double h2)
{  return (1.0/6.0)*h2*(ib3(x, klo, khi, tau, h, h2) - ib(x, klo, khi, tau, h, h2)); }



double mlfi_spline_fi(const spline *s, const double xi)
{
  unsigned long klo, khi;
  double h, h2, iad, ibd, icd, idd;

  klo = find(s, xi);
  khi = klo + 1L;

  h = s->x[khi] - s->x[klo];
  h2 = h*h;
  iad = ia(s->x, klo, khi, xi, h, h2) - ia(s->x, klo, khi, s->x[klo], h, h2);
  ibd = ib(s->x, klo, khi, xi, h, h2) - ib(s->x, klo, khi, s->x[klo], h, h2);
  icd = ic(s->x, klo, khi, xi, h, h2) - ic(s->x, klo, khi, s->x[klo], h, h2);
  idd = id(s->x, klo, khi, xi, h, h2) - id(s->x, klo, khi, s->x[klo], h, h2);

  return s->fc[klo] + iad*s->y[klo] + ibd*s->y[khi] + icd*s->yxx[klo] + idd*s->yxx[khi];
}

void tridiag_lu_init(int n, double *al, double *ml)
{
  double *ad = al + n;
  double *au = ad + n;
  double *md = ml + n;
  double *mu = md + n;
  int j;

  ad[0] = md[0];
  for (j=1; j < n; j++) {
    au[j-1] = mu[j-1];
    al[j]   = ml[j] / ad[j-1];
    ad[j]   = md[j] - al[j] * au[j-1];
  }
}

void tridiag_lu_solve(int n,
                      double *al,
                      double *w, double *v, double *q)
{
  double *ad = al + n;
  double *au = ad + n;
  int j;

  w[0] = q[0];
  for (j=1; j < n; j++) {
    w[j] = q[j] - al[j] * w[j-1];
  }

  v[n-1] = w[n-1] / ad[n-1];
  for (j=n-2; 0 <= j; j--) {
    v[j] = (w[j] - au[j] * v[j+1]) / ad[j];
  }
}

void tridiag_lu_init_solve (const int n, double *a, const double *v, double *x)
{
  double *b = a + n;
  double *c = b + n;
  int i;

  x[0] = v[0];
  for (i = 1; i < n; i++)
    {
      double m = a[i] / b[i-1];
      b[i] = b[i] - m * c[i-1];
      x[i] = v[i] - m * x[i-1];
    }

  x[n-1] = x[n-1] / b[n-1];
  for (i = n - 2; i >= 0; i--)
    x[i] = (x[i] - c[i] * x[i+1]) / b[i];
}

void interior_scheme_1d (const int n, double *ml, const double *a_xx, const double *a_x, const double *a, const double dx, const double dt)
{
  double *md = ml + n;
  double *mu = md + n;
  double coeff_a_x = dt / dx;
  double coeff_a_xx = coeff_a_x / dx;
  int i;
  /* left = dt / dx^2 * a_xx - 0.5 * dt / dx * a_x */
  /* center = 1. + dt * a - 2.* dt / dx^2 * a_xx */
  /* right = dt / dx^2 * a_xx + 0.5 * dt / dx * a_x */
  for (i = 1; i < n - 1; i++)
    {
      double mod_axxi = a_xx[i] * coeff_a_xx;
      double mod_axi = a_x[i] * coeff_a_x;
      double mod_ai = a[i] * dt;
      ml[i] = mod_axxi - 0.5 * mod_axi;
      md[i] = 1.0 + mod_ai - 2.0 * mod_axxi;
      mu[i] = mod_axxi + 0.5 * mod_axi;
    }
}

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
                                      const int n)
{
  int i;
  for (i = 0; i < n; i++) {
    s[i] = ft * gx[i] + kt;
    j_x[i] = ft * dgx[i];
    h_x[i] = ft * d2gx[i];
    j_t[i] = dft * gx[i] + dkt;
  }
}

void apply_change_of_variable_1d (double* a_xx, double* a_x, const double* j_x, const double* j_t, const double* h_x , const int n)
{
  int i;
  double j, mod_axx;
  for (i = 0; i < n; i++) {
    j = j_x[i];
    mod_axx = a_xx[i] / (j * j);
    a_xx[i] = mod_axx;
    a_x[i] = (a_x[i] - j_t[i] - mod_axx * h_x[i]) / j;
  }
}

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
                                  const int ny)
{
  int i, j, idx;
  double jac, hess, jact, mod_b, coeff_b, inv_j;
  for (i = 0; i < nx; i++) {
    jac = j_x[i];
    hess = h_x[i];
    jact = jx_t[i];
    inv_j = 1 / jac;
    coeff_b = inv_j * inv_j;
    for (j = 0; j < ny; j++) {
      idx = i * ny + j;
      mod_b = coeff_b * b_xx[idx];
      b_xx[idx] = mod_b;
      b_xy[idx] *= inv_j / j_y[j];
      b_x[idx] = (b_x[idx] - jact - mod_b * hess) * inv_j;
    }
  }
  for (j = 0; j < ny; j++) {
    jac = j_y[j];
    hess = h_y[j];
    jact = jy_t[j];
    inv_j = 1 / jac;
    coeff_b = inv_j * inv_j;
    for (i = 0; i < nx; i++) {
      idx = i * ny + j;
      mod_b = coeff_b * b_yy[idx];
      b_yy[idx] = mod_b;
      b_y[idx] = (b_y[idx] - jact - mod_b * hess) * inv_j;
    }
  }
}

void tridiag_gemv(enum CBLAS_TRANSPOSE v_trans,
                  const int n,
                  double *al,
                  const double alpha,
                  double *x,
                  const double beta,
                  double *y)
{
  double *ad = al + n;
  double *au = ad + n;
  double *temp;
  int i;

  if (v_trans == CblasTrans) {temp = au; au = al+1; al = temp-1;};
  /*
    al = xx l2 l3 l4 l5 l6
    ad = d1 d2 d3 d4 d5 d6
    au = u1 u2 u3 u4 u5 yy
   */
  y[0] = alpha * (ad[0] * x[0] + au[0] * x[1]) + beta * y[0];
  for (i = 1; i < n-1; i++) {
    y[i] = alpha * (al[i] * x[i-1] + ad[i] * x[i] + au[i] * x[i+1]) + beta * y[i];
  };
  y[n-1] = alpha * (al[n-1] * x[n-2] + ad[n-1] * x[n-1]) + beta * y[n-1];
}
