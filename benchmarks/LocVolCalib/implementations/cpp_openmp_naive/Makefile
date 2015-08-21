# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include

HELPERS     =$(HIPERMARK_BENCHMARK_LIB_DIR)/include/ParseInput.h
OBJECTS     =VolCalibOrig.o
EXECUTABLE  =VolCalibOrig

default: $(EXECUTABLE)

%.o: $(HIPERMARK_IMPLEMENTATION_DIR)/%.cpp $(HELPERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

data: $(HIPERMARK_INPUT)
	mkdir -p datasets/$(HIPERMARK_INPUT_NAME)
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT) outer num_x num_y num_t > datasets/$(HIPERMARK_INPUT_NAME)/input.data

run:
	cat datasets/$(HIPERMARK_INPUT_NAME)/input.data | \
	 OMP_NUM_THREADS=$(HIPERMARK_CONFIG_NUM_THREADS) ./$(EXECUTABLE)
