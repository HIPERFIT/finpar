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

HIPERMARK_DATA_FIELDS=outer num_x num_y num_t
HIPERMARK_RUN_ENVIRONMENT=
include $(HIPERMARK_LIB_DIR)/flat_data_helper.mk
