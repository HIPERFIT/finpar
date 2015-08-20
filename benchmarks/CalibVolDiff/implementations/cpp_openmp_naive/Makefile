# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

include $(HIPERMARK_LIB_DIR)/setup.mk

SOURCES_CPP =VolCalibOrig.cpp
HELPERS     =$(HIPERMARK_BENCHMARK_LIB_DIR)/ParseInput.h
OBJECTS     =VolCalibOrig.o
EXECUTABLE  =VolCalibOrig

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

cpu: $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

data: $(HIPERMARK_INPUT)
	mkdir -p datasets/$(HIPERMARK_INPUT_NAME)
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT) outer num_x num_y num_t > datasets/$(HIPERMARK_INPUT_NAME)/input.data
