# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

# FIXME - this should not even exist!
include $(HIPERMARK_LIB_DIR)/platform.mk

# Handle Strength reduction and vectorization
ifeq ($(HIPERMARK_CONFIG_STRENGTH_REDUCTION),on)
	CXXFLAGS += -D_OPTIMIZATION_SOBOL_STRENGTH_RED_RECURR
endif
ifeq ($(HIPERMARK_CONFIG_VECTORIZATION),auto)
	CXXFLAGS += -D_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV=0
else ifeq ($(HIPERMARK_CONFIG_VECTORIZATION),off)
	CXXFLAGS += -D_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV=2
else ifeq ($(HIPERMARK_CONFIG_VECTORIZATION),on)
	CXXFLAGS += -D_OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV=1
endif



INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include

GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) -D GPU_DEV_ID=$(GPU_DEVICE_ID) \
             -D GPU_LOC_MEM=$(GPU_LOCAL_MEM) -D GPU_CONST_MEM=$(GPU_CONST_MEM) \
             -D GPU_REG_MEM=$(GPU_REG_MEM) -D GPU_GLB_MEM=$(GPU_DEVICE_MEM) \
             -D GPU_TILE=$(GPU_LOCAL_MEM_PER_TH) -D GPU_CORES=$(GPU_NUM_CORES) \
		-D CURR_DIR_PATH='"$(MAKE_DIR)"'

SOURCES_CPP =GenPricing.cpp
HELPERS     =$(HIPERMARK_BENCHMARK_LIB_DIR)/include/ParseInput.h
OBJECTS     =GenPricing.o
EXECUTABLE  =GenPricing

default: $(EXECUTABLE)

%.o: $(HIPERMARK_IMPLEMENTATION_DIR)/%.cpp $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)
	@# clean nVidia compiler cache
	rm -rf $(HOME)/.nv/ComputeCache/*

HIPERMARK_DATA_FIELDS=contract_number monte_carlo_iterations num_path_dates num_underlyings num_models num_bits dir_vs md_c md_vols md_drifts md_starts md_deterministic_values md_discounts bb_inds bb_data
HIPERMARK_RUN_ENVIRONMENT=
include $(HIPERMARK_LIB_DIR)/flat_data_helper.mk
