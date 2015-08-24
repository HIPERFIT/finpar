# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

include $(HIPERMARK_LIB_DIR)/setup.mk

SOURCES    = GenPricing.hs
EXECUTABLE = GenPricing
OBJECTS    = GenPricing.hi GenPricing.o

GHC_FLAGS  = -Wall -O2 -msse2 -rtsopts -i$(HIPERMARK_IMPLEMENTATION)

default: $(EXECUTABLE)

$(EXECUTABLE): $(HIPERMARK_IMPLEMENTATION_DIR)/$(SOURCES)
	ghc $(HIPERMARK_IMPLEMENTATION_DIR)/$(SOURCES) $(GHC_FLAGS) -o $(EXECUTABLE)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

data: $(HIPERMARK_INPUT)
	mkdir -p datasets/$(HIPERMARK_INPUT_NAME)
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT) contract_number monte_carlo_iterations num_path_dates num_underlyings num_models num_bits dir_vs md_c md_vols md_drifts md_starts md_deterministic_values md_discounts bb_inds bb_data > datasets/$(HIPERMARK_INPUT_NAME)/input.data

run:
	cat datasets/$(HIPERMARK_INPUT_NAME)/input.data | \
	 ./$(EXECUTABLE)
