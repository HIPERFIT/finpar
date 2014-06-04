/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: mt19937ar.c 41326 2012-01-03 17:30:52Z jmeber $ */

/* 
   Mersenne Twister random number generator. With small modifications
   by LexiFi SAS to replace global state with a structure. 
   
*/


/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

#include "mt19937ar.h"

/* Period parameters */  
#define N MERSENNE_CYCLE
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* initializes mt[N] with a seed */
void mt_init_genrand(rng *mtp, unsigned long s)
{
    mtp->mt[0]= s & 0xffffffffUL;
    for (mtp->mti=1; mtp->mti<N; mtp->mti++) {
        mtp->mt[mtp->mti] = 
	    (1812433253UL * (mtp->mt[mtp->mti-1] ^ (mtp->mt[mtp->mti-1] >> 30)) + mtp->mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mtp->mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mtp->mt[mtp->mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void mt_init_by_array(rng *mtp, unsigned long init_key[],int key_length)
{
    int i, j, k;
    mt_init_genrand(mtp, 19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mtp->mt[i] = (mtp->mt[i] ^ ((mtp->mt[i-1] ^ (mtp->mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mtp->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mtp->mt[0] = mtp->mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mtp->mt[i] = (mtp->mt[i] ^ ((mtp->mt[i-1] ^ (mtp->mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mtp->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mtp->mt[0] = mtp->mt[N-1]; i=1; }
    }

    mtp->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_genrand_int32(rng *mtp)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mtp->mti >= N) { /* generate N words at one time */
        int kk;

        if (mtp->mti == N+1)   /* if init_genrand() has not been called, */
            mt_init_genrand(mtp, 5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mtp->mt[kk]&UPPER_MASK)|(mtp->mt[kk+1]&LOWER_MASK);
            mtp->mt[kk] = mtp->mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mtp->mt[kk]&UPPER_MASK)|(mtp->mt[kk+1]&LOWER_MASK);
            mtp->mt[kk] = mtp->mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mtp->mt[N-1]&UPPER_MASK)|(mtp->mt[0]&LOWER_MASK);
        mtp->mt[N-1] = mtp->mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mtp->mti = 0;
    }
  
    y = mtp->mt[mtp->mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long mt_genrand_int31(rng *mtp)
{
  return (long)(mt_genrand_int32(mtp)>>1);
}

/* generates a random number on [0,1]-real-interval */
double mt_genrand_real1(rng *mtp)
{
  return mt_genrand_int32(mtp)*(1.0/4294967295.0); 
  /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double mt_genrand_real2(rng *mtp)
{
    return mt_genrand_int32(mtp)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double mt_genrand_real3(rng *mtp)
{
    return (((double)mt_genrand_int32(mtp)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double mt_genrand_res53(rng *mtp) 
{ 
    unsigned long a=mt_genrand_int32(mtp)>>5, b=mt_genrand_int32(mtp)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

int testmain(void)
{
    int i;
    rng mtp; 
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    
    mtp.mti=N+1;

    mt_init_by_array(&mtp, init, length);
    printf("1000 outputs of genrand_int32()\n");
    for (i=0; i<1000; i++) {
      printf("%10lu ", mt_genrand_int32(&mtp));
      if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand_real2()\n");
    for (i=0; i<1000; i++) {
      printf("%10.8f ", mt_genrand_real2(&mtp));
      if (i%5==4) printf("\n");
    }
    return 0;
}
