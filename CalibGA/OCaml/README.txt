###########################################################################
#  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               #
#                                                                         #
#  No part of this document may be reproduced or transmitted in any       #
#  form or for any purpose without the express permission of LexiFi SAS.  #
###########################################################################

# $Id: README.txt 1158 2012-12-06 13:59:49Z afrisch $

This directory contains a simplified version of LexiFi's calibration
routine for the G2++ model.

Usage: simply type 'make'.  This will compile and run a test program.
The code requires OCaml 4.00 or above.

The calibration routine is based on least-square optimization, using
the DE algorithm (in LexiFi's version, this global optimizer is
followed by a local one).  The optimizer is implemented in the
Optimization module.

The Math module provides:
 - an approximation formula for the cumulative distribution function
   for a normal variate
 - a root finder (using Brent's algorithm)
 - coefficients for a Gaussian quadrature (for quick integral based
   on sampling a few points of the integrand)

The G2pp_calibration module implements to core calibration routine.
