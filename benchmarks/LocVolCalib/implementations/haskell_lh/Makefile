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

data: $(HIPERMARK_INPUT)
	mkdir -p datasets/$(HIPERMARK_INPUT_NAME)
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT) outer num_x num_y num_t > datasets/$(HIPERMARK_INPUT_NAME)/input.data

run:
	cat datasets/$(HIPERMARK_INPUT_NAME)/input.data | ./$(EXECUTABLE)
