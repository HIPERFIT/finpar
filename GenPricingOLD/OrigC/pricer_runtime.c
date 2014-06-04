/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* $Id: pricer_runtime.c 41326 2012-01-03 17:30:52Z jmeber $ */

/*
  MLFi pricer code runtime library.
  These routines must be made available to any MLFi generated pricer.
  A generated pricer may call any of these routines.
 */

#include "mlfi_pricing.h"
#include "mlfi_erfc.h"

/* general float routines */


double erfc(const double x)
{return mlfi_erfc(x);}

/* array ordering, sorting, extraction. */

long array_index(const int n, const int elm, const int* arr)
{
  int res, i;
  res = -1;
  for (i=0;i<n;i++){
    if (arr[i] == elm) res = i;
  };
  return res;
}

void order(int* dest, const int n, const int* arr, const double* fs)
{
  int i,j;
  int a;
  double fa;
  dest[0]=arr[0];
  for (j=1;j<n;j++){
    a=arr[j];
    fa=fs[a];
    i=j-1;
    while (i >= 0 && fs[dest[i]] > fa){
      dest[i+1]=dest[i];
      i--;}
    dest[i+1]=a;}
}
