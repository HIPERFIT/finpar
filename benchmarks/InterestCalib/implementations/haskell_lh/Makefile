# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

include $(HIPERMARK_LIB_DIR)/setup.mk

SOURCES    = SwapCalib.hs
EXECUTABLE = SwapCalib
OBJECTS    = SwapCalib.hi SwapCalib.o

GHC_FLAGS  = -Wall -O2 -msse2 -rtsopts -i$(HIPERMARK_IMPLEMENTATION_DIR)

default: $(EXECUTABLE)

$(EXECUTABLE): $(HIPERMARK_IMPLEMENTATION_DIR)/$(SOURCES)
	ghc $(HIPERMARK_IMPLEMENTATION_DIR)/$(SOURCES) $(GHC_FLAGS) -o $(EXECUTABLE)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

data: $(HIPERMARK_INPUT)
	mkdir -p datasets/$(HIPERMARK_INPUT_NAME)
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT) population_size convergence_iterations num_swaption_quotes swaption_quotes num_hermite_coefficients hermite_coefficients hermite_coefficient_weights sobol_num_bits sobol_dir_vs > datasets/$(HIPERMARK_INPUT_NAME)/input.data

run:
	cat datasets/$(HIPERMARK_INPUT_NAME)/input.data | \
	 ./$(EXECUTABLE)
