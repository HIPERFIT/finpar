# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

# FIXME - this should not even exist!
include $(HIPERMARK_LIB_DIR)/platform.mk

INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include

GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) -D GPU_DEV_ID=$(GPU_DEVICE_ID) \
             -D GPU_LOC_MEM=$(GPU_LOCAL_MEM) -D GPU_CONST_MEM=$(GPU_CONST_MEM) \
             -D GPU_REG_MEM=$(GPU_REG_MEM) -D GPU_GLB_MEM=$(GPU_DEVICE_MEM) \
             -D GPU_TILE=$(GPU_LOCAL_MEM_PER_TH) -D GPU_CORES=$(GPU_NUM_CORES) \
			 -D CURR_DIR_PATH='"$(MAKE_DIR)"'

SOURCES_CPP =SwapCalib.cpp
HELPERS     =$(HIPERMARK_BENCHMARK_LIB_DIR)/include/ParseInput.h
OBJECTS     =SwapCalib.o
EXECUTABLE  =SwapCalib


default: $(EXECUTABLE)

%.o: $(HIPERMARK_IMPLEMENTATION_DIR)/%.cpp $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

$(EXECUTABLE): $(OBJECTS)
	cp -r $(HIPERMARK_IMPLEMENTATION_DIR)/SrcCL .
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

run:
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT_FILE) population_size convergence_iterations num_swaption_quotes swaption_quotes num_hermite_coefficients hermite_coefficients hermite_coefficient_weights sobol_num_bits sobol_dir_vs | \
	 ./$(EXECUTABLE)
