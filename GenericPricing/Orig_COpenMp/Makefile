include ../../setup.mk

INCLUDES    += -I../includeC

SOURCES_CPP =GenPricing.cpp
HELPERS     =Contracts.h SobolGaussBB.h
# XXX what about the boilerplate in ../includeC ???
OBJECTS     =GenPricing.o
EXECUTABLE  =GenPricing

default: cpu

# TODO need redesign, probably use different executable names (or
# standardisation of "gpu" targets across all benchmarks)

#.cpp.o:
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

gpu:
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

cpu: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) -D IS_GPU=0 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP)

default: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_medium: $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_large: $(EXECUTABLE)
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE) 2> Debug.txt

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

