# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_IMPLEMENTATION_DIR)

#-fno-unsafe-math-optimizations -fno-finite-math-only 
#-fmath-errno -ftrapping-math -fno-rounding-math -fsignaling-nans

HELPERS     =$(HIPERMARK_BENCHMARK_LIB_DIR)/include/ParseInput.h
OBJECTS     =SwapCalib.o
EXECUTABLE  =SwapCalib

default: $(EXECUTABLE)

%.o: $(HIPERMARK_IMPLEMENTATION_DIR)/%.cpp $(HELPERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

run:
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT_FILE) population_size convergence_iterations num_swaption_quotes swaption_quotes num_hermite_coefficients hermite_coefficients hermite_coefficient_weights sobol_num_bits sobol_dir_vs | \
	 OMP_NUM_THREADS=$(HIPERMARK_CONFIG_NUM_THREADS) ./$(EXECUTABLE)
