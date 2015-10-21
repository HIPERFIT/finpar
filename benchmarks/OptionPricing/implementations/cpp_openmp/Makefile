# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include

SOURCES_CPP =GenPricing.cpp
HELPERS     =$(HIPERMARK_BENCHMARK_LIB_DIR)/include/ParseInput.h
OBJECTS     =GenPricing.o
EXECUTABLE  =GenPricing

default: $(EXECUTABLE)

%.o: $(HIPERMARK_IMPLEMENTATION_DIR)/%.cpp $(HELPERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

HIPERMARK_DATA_FIELDS=contract_number monte_carlo_iterations num_path_dates num_underlyings num_models num_bits dir_vs md_c md_vols md_drifts md_starts md_deterministic_values md_discounts bb_inds bb_data
ifdef HIPERMARK_CONFIG_NUM_THREADS
  HIPERMARK_RUN_ENVIRONMENT=OMP_NUM_THREADS=$(HIPERMARK_CONFIG_NUM_THREADS)
endif
include $(HIPERMARK_LIB_DIR)/flat_data_helper.mk
