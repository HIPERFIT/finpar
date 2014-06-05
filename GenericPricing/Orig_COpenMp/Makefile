OPENCL_ROOTDIR    ?= /usr/local/cuda
OPENCL_LIBDIR     := $(OPENCL_ROOTDIR)/lib64
OPENCL_INCDIR	  ?= $(OPENCL_ROOTDIR)/include

CXX        = g++
#INCLUDES   = -I$(CUDA_INCDIR) -I. -I../includeC -I../../include
INCLUDES   = -I$(OPENCL_INCDIR) -I. -I../includeC -I../../include
# TODO CUDA or OPENCL here?
CXXFLAGS   = -O2 -fno-strict-aliasing -fopenmp
OPENCL     = -L$(OPENCL_LIBDIR) -L/usr/lib/ -lOpenCL
# TODO is /usr/lib really needed here?
#-m64
#-ffast-math  -fopenmp

SOURCES_CPP =GenPricing.cpp
HELPERS     =Contracts.h SobolGaussBB.h
# XXX what about the boilerplate in ../includeC ???
OBJECTS     =GenPricing.o
EXECUTABLE  =GenPricing

#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL

default: cpu

# TODO need redesign, probably use different executable names (or
# standardisation of "gpu" targets across all benchmarks)

#.cpp.o:
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

gpu:
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(OPENCL)

cpu: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) -D IS_GPU=0 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP)

run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE)


run_medium:  $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE)

run_large:  $(EXECUTABLE)
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

