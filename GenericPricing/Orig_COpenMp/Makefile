OPENCL_ROOTDIR    ?= /usr/local/cuda
OPENCL_LIBDIR     := $(OPENCL_ROOTDIR)/lib64
OPENCL_INCDIR	  ?= $(OPENCL_ROOTDIR)/include

CXX        = g++
INCLUDES   = -I$(CUDA_INCDIR) -I. -I../includeC -I../../include
LIB        = -L$(OPENCL_LIBDIR) -L/usr/lib/ -lOpenCL
CXXFLAGS   = -O2 -fno-strict-aliasing -fopenmp
#-m64
#-ffast-math  -fopenmp

SOURCES_CPP =GenPricing.cpp
OBJECTS     =GenPricing.o  
EXECUTABLE  =GenPricing

#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL


.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

gpu: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

cpu: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=0 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

default: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

run_small: 
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE)

run_medium: 
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE)

run_large: 
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE)

#default: $(OBJECTS)
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

