OPENCL_ROOTDIR    ?= /usr/local/cuda
OPENCL_LIBDIR     := $(OPENCL_ROOTDIR)/lib64
OPENCL_INCDIR	  ?= $(OPENCL_ROOTDIR)/include

CXX        = g++
INCLUDES   = -I$(CUDA_INCDIR) -I.
LIB        = -L$(OPENCL_LIBDIR) -L/usr/lib/ -lOpenCL
CXXFLAGS   = -O3 -fopenmp

SOURCES_CPP =NordeaVect.cpp
OBJECTS     =NordeaVect.o  
EXECUTABLE  =nordea

#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL


.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

gpu: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

cpu: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=0 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

default: 
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

#default: $(OBJECTS)
#	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
