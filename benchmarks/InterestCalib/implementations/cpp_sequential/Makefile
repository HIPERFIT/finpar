# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include

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

HIPERMARK_DATA_FIELDS=population_size convergence_iterations num_swaption_quotes swaption_quotes num_hermite_coefficients hermite_coefficients hermite_coefficient_weights sobol_num_bits sobol_dir_vs
HIPERMARK_RUN_ENVIRONMENT=
include $(HIPERMARK_LIB_DIR)/flat_data_helper.mk
