include ../../setup.mk

ifndef ($(HAVE_GPU))
include ../../platform.mk
endif

INCLUDES    += -I../includeC
GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) 

SOURCES_CPP =GenPricing.cpp
HELPERS     =Contracts.h SobolGaussBB.h ../includeC/Constants.h ../includeC/Optimizations.h ../includeC/ParseInput.h
# XXX what about the boilerplate in ../includeC ???
OBJECTS     =GenPricing.o
EXECUTABLE  =GenPricing

default: cpu

# TODO need redesign, probably use different executable names (or
# standardisation of "gpu" targets across all benchmarks)

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

cpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_medium: $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_large: $(EXECUTABLE)
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE) 2> Debug.txt

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

