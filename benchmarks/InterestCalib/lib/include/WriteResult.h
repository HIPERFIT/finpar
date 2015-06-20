#ifndef WRITERESULT_H
#define WRITERESULT_H

#include <stdio.h>

/* Writes trailing comma and newline. */
void writeRealField(FILE *stream, const char* field, REAL x)
{
  fprintf(stream, "\"%s\": %f,\n", field, x);
}

void writeResult(FILE *stream,
                 REAL a,
                 REAL b,
                 REAL sigma,
                 REAL nu,
                 REAL rho,
                 REAL logLik,
                 REAL* calib_arr)
{
  /* Hackily manually generate JSON.  It would be nice with some
     proper library support. */
  fprintf(stream, "{\n");

  writeRealField(stream, "a_field", a);
  writeRealField(stream, "b_field", b);
  writeRealField(stream, "sigma_field", sigma);
  writeRealField(stream, "nu_field", nu);
  writeRealField(stream, "rho_field", rho);
  writeRealField(stream, "lg_likelyhood", logLik);

  fprintf(stream, "\"swaption_calibration_result\":\n");
  write_2Darr(stream, calib_arr, static_cast<int>(NUM_SWAP_QUOTES), 3);
  fprintf(stream, "}\n");
}

#endif
