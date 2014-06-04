/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mlfi_sobol.h 41326 2012-01-03 17:30:52Z jmeber $ */

#include "mt19937ar.h"

typedef struct
{
  const char * name;
  unsigned int max_dimension;
  size_t (*state_size) (unsigned int dimension);
  void (*init_state) (void * state, unsigned int dimension);
  void (*get) (void * state, unsigned int dimension, double x[]);
} qrng_type;


typedef struct
{
  const qrng_type * type;
  unsigned int dimension;
  size_t state_size;
  void * state;
} qrng;

qrng_type *qrng_sobol();

qrng *qrng_alloc (const qrng_type * T, unsigned int dimension);
void qrng_free (qrng *r);
void qrng_get (const qrng * r, double x[]);

#define SOBOL_BIT_COUNT 30

/* Sobol generator state.
 *   sequence_count       = number of calls with this generator
 *   last_numerator_vec   = last generated numerator vector
 *   last_denominator_inv = 1/denominator for last numerator vector
 *   v_direction          = direction number table
 */
typedef struct
{
  unsigned int  sequence_count;
  double        last_denominator_inv;
  int           *last_numerator_vec;
  int           *v_direction[SOBOL_BIT_COUNT];
  rng           r;
} sobol_state_t;

sobol_state_t *sobol_create(int dimension);
void sobol_init(void *state, unsigned int dimension);
void sobol_get(void *state, unsigned int dimension, double * v);
void sobol_free(sobol_state_t *state);
