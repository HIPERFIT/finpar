/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mlfi_sobol.c 42155 2012-01-27 13:52:19Z jmeber $ */

#include <stdio.h>
#include <stdlib.h>
#include "mt19937ar.h"
#include "mlfi_sobol.h"
#include "poly.h"

/* prototypes for generator type functions */
static size_t sobol_state_size(unsigned int dimension);

/* global Sobol generator type object */
static const qrng_type sobol_type =
{
  "sobol",
  0,
  sobol_state_size,
  sobol_init,
  sobol_get
};

qrng_type * qrng_sobol()
{
  return (qrng_type*) &sobol_type;
}

static size_t sobol_state_size(unsigned int dimension)
{
  return
    sizeof(sobol_state_t) +      /* The struct */
    sizeof(int) * dimension +    /* for last_numerator_vec */
    sizeof(int) * dimension * SOBOL_BIT_COUNT; /* for the direction no.s */
}

unsigned long int uniform_int (rng * r, unsigned long int n)
{
  unsigned long int offset = 0UL;
  unsigned long int range = 0xffffffffUL  - offset;
  unsigned long int scale = range / n;
  unsigned long int k;

  do {
    k = (mt_genrand_int32(r) - offset) / scale;
  }
  while (k >= n);

  return k;
}

/* l is a degree number, between 1 and the degree of the polynomial
   associated with the dimension being initialized minus 1 */
static int sobol_init_direction(rng *r, int l)
{
  /* See Peter Jaeckel, Monte Carlo Methods in Finance, Wiley 2002, p86.
   */
  int wkl = uniform_int(r, 1<<l) | 1;
  /*  printf("l=%d n=%ld wkl=%d mask=%ld\n", l, n, wkl, mask);*/
  return wkl;
}

static void failure(char *str)
{
  fprintf(stderr, "sobol: %s\n", str);
  exit(1);
}

void sobol_init(void * state, unsigned int dimension)
{
  sobol_state_t * s_state = (sobol_state_t *) state;
  unsigned int i_dim, dim;
  int j, k;
  int ell;
  int *includ;
  int max_degree;

  if(dimension < 0)
    failure("dimension should be >= 0");

  dim = dimension > max_poly ? max_poly : dimension;

  if (dimension > max_poly) {
    fprintf(stderr, "warning: sobol dimension limit exceeded, using mersenne rng\n");
  }

  mt_init_genrand(&s_state->r, 0L);

  if (dim > 0) { /* initialize includ */

    max_degree = degrees[dim-1];
    includ = (int *) malloc(sizeof(int) * max_degree);
    if (includ==0) {
      failure ("allocation of 'includ' failed for sobol init");
    }

    /* Initialize direction table in dimension 0. */
    for(k=0; k<SOBOL_BIT_COUNT; k++) {
      s_state->v_direction[k] =
	(int *) (((char *) s_state) +
		 sizeof(sobol_state_t) +
		 sizeof(int) * dim +
		 sizeof(int) * dim * k);
      s_state->v_direction[k][0] = 1;
    }

  }
  else {
    max_degree=0;
    includ=0;
  }

  s_state->last_numerator_vec = (int *) ((char *) s_state + sizeof(sobol_state_t));

  /* Initialize in remaining dimensions. */
  for(i_dim=1; i_dim<dim; i_dim++) {

    const int poly_index = i_dim;
    const int degree_i = degrees[poly_index];

    /* Expand the polynomial bit pattern to separate
     * components of the logical array includ[].
     */
    int p_i = polynomials[poly_index];
    for(k = degree_i-1; k >= 0; k--) {
      includ[k] = ((p_i % 2) == 1);
      p_i /= 2;
    }

    /* Leading elements for dimension i are randomly initialized */
    for(j=0; j<degree_i; j++) s_state->v_direction[j][i_dim] = sobol_init_direction(&s_state->r, j+1);

    /* Calculate remaining elements for this dimension,
     * as explained in Bratley+Fox, section 2.
     */
    for(j=degree_i; j<SOBOL_BIT_COUNT; j++) {
      int newv = s_state->v_direction[j-degree_i][i_dim];
      ell = 1;
      for(k=0; k<degree_i; k++) {
        ell *= 2;
        if(includ[k]) newv ^= (ell * s_state->v_direction[j-k-1][i_dim]);
      }
      s_state->v_direction[j][i_dim] = newv;
    }
  }

  /* Multiply columns of v by appropriate power of 2. */
  ell = 1;
  for(j=SOBOL_BIT_COUNT-1-1; j>=0; j--) {
    ell *= 2;
    for(i_dim=0; i_dim<dim; i_dim++) {
      s_state->v_direction[j][i_dim] *= ell;
    }
  }


  /* 1/(common denominator of the elements in v_direction) */
  s_state->last_denominator_inv = 1.0 /(2.0 * ell);

  /* final setup */
  s_state->sequence_count = 0;
  for(i_dim=0; i_dim<dim; i_dim++) s_state->last_numerator_vec[i_dim] = 0;

  if (dim > 0)
    free(includ);
}


#define SOBOL_BIT_TABLE 1

void sobol_get(void * state, unsigned int dimension, double * v)
{
	int jj;
  sobol_state_t * s_state = (sobol_state_t *) state;
  unsigned int i_dimension, dim;

  /* Find the position of the least-significant zero in count. */

#ifdef SOBOL_BIT_TABLE

#define L1 1
#define L2 2, L1
#define L3 3, L1, L2
#define L4 4, L1, L2, L3
#define L5 5, L1, L2, L3, L4
#define L6 6, L1, L2, L3, L4, L5
#define L7 7, L1, L2, L3, L4, L5, L6
#define L8 8, L1, L2, L3, L4, L5, L6, L7

  register unsigned int c = ~(s_state->sequence_count);
  int ell;
  register unsigned int d;

  static const int BitTable[256] =
    {
      0, L1, L2, L3, L4, L5, L6, L7, L8
    };
#if 0
  printf("Printing Bit Table: \n\t");
  for(jj=0; jj<256; jj++) {
	  printf("%d, ", BitTable[jj]);
  }
  printf("\n\n");
#endif

  if (d = c & 0xff)
    ell = BitTable[d];
  else if (d = (c >> 8) & 0xff)
    ell = BitTable[d] + 8;
  else if (d = (c >> 16) & 0xff)
    ell = BitTable[d] + 16;
  else
    ell = BitTable[c >> 24] + 24;

#else

  int ell = 1;
  unsigned int c = s_state->sequence_count;

  while(c & 1) {
    ell++;
    c >>= 1;
  }
#endif

  /* Check for exhaustion. */
  if(ell > SOBOL_BIT_COUNT) failure("bit count exceeded");

  dim = dimension > max_poly ? max_poly : dimension;

  for(i_dimension=0; i_dimension<dim; i_dimension++) {
    const int direction_i     = s_state->v_direction[ell-1][i_dimension];
    const int old_numerator_i = s_state->last_numerator_vec[i_dimension];
    const int new_numerator_i = old_numerator_i ^ direction_i;
    s_state->last_numerator_vec[i_dimension] = new_numerator_i;
    v[i_dimension] = new_numerator_i * s_state->last_denominator_inv;
  }
  for (i_dimension = dim; i_dimension < dimension; i_dimension++) {
    v[i_dimension] = mt_genrand_real2(&s_state->r);
  }

  s_state->sequence_count++;

  return;
}

void sobol_free(sobol_state_t *state) {
  free(state);
}

sobol_state_t *sobol_create(int dimension) {
  sobol_state_t *s;
  s = malloc(sobol_state_size(dimension));
  sobol_init(s, dimension);
  return s;
}

qrng *qrng_alloc (const qrng_type * T, unsigned int dimension)
{

  qrng *r = (qrng *) malloc (sizeof (qrng));

  if (r == 0) failure("allocation failed for qrng struct");

  r->dimension = dimension;
  r->state_size = T->state_size(dimension);
  r->state = malloc (r->state_size);

  if (r->state == 0) {
    free (r);
    failure("allocation failed for qrng state");
  }

  r->type = T;

  T->init_state(r->state, r->dimension);

  return r;
}

void qrng_free (qrng *r)
{
  if (r != 0) {
    if (r->state != 0) free(r->state);
    free(r);
  }
}

void qrng_get (const qrng * r, double x[])
{
  (r->type->get) (r->state, r->dimension, x);
}
