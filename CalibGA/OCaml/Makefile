###########################################################################
#  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               #
#                                                                         #
#  No part of this document may be reproduced or transmitted in any       #
#  form or for any purpose without the express permission of LexiFi SAS.  #
###########################################################################

# $Id: Makefile 1158 2012-12-06 13:59:49Z afrisch $

FILES=\
  date.mli \
  date.ml \
  math.mli \
  math.ml \
  optimization.mli \
  optimization.ml \
  g2pp_calibration.mli \
  g2pp_calibration.ml \
  main.ml

run:
	ocamlopt -o run.exe -w +A $(FILES)
	time ./run.exe

pack:
	tar czf g2pp.tar.gz $(FILES) Makefile README.txt

clean:
	rm -rf *.cm* run.exe *.o *.obj *~ *.tar.gz
