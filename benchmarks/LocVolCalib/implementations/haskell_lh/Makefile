# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

include $(HIPERMARK_LIB_DIR)/setup.mk

SOURCES    = VolCalib.hs
EXECUTABLE = VolCalib
GHC_FLAGS  = -Wall -O2 -msse2 -rtsopts

default: $(EXECUTABLE)

$(EXECUTABLE): $(HIPERMARK_IMPLEMENTATION_DIR)/$(SOURCES)
	ghc $(HIPERMARK_IMPLEMENTATION_DIR)/$(SOURCES) $(GHC_FLAGS) -o $(EXECUTABLE)

clean:
	rm $(EXECUTABLE) *o *hi

HIPERMARK_DATA_FIELDS=outer num_x num_y num_t
HIPERMARK_RUN_ENVIRONMENT=
include $(HIPERMARK_LIB_DIR)/flat_data_helper.mk
