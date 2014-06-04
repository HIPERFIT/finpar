/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

#ifndef _MT19937AR_H
#define _MT19937AR_H

#define MERSENNE_CYCLE 624

typedef struct {
  unsigned long mt[MERSENNE_CYCLE]; /* the array for the state vector  */
  int mti;          /* mti==MERSENNE_CYCLE+1 means mt[MERSENNE_CYCLE] is not initialized */
} rng;

void mt_init_genrand(rng *mtp, unsigned long s);
void mt_init_by_array(rng *mtp, unsigned long init_key[],int key_length);

double mt_genrand_real1(rng *mtp); /* on 0 <= x <= 1 */
double mt_genrand_real2(rng *mtp); /* on 0 <= x <  1 */
double mt_genrand_real3(rng *mtp); /* on 0 <  x <  1 */
double mt_genrand_res53(rng *mtp); /* on 0 <= x <  1 with 53 bits */

unsigned long mt_genrand_int32(rng *mtp);
long mt_genrand_int31(rng *mtp);

#endif
