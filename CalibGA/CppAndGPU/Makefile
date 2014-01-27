OPENCL_ROOTDIR    ?= /usr/local/cuda
OPENCL_LIBDIR     := $(OPENCL_ROOTDIR)/lib64
OPENCL_INCDIR	  ?= $(OPENCL_ROOTDIR)/include

CXX        = g++
INCLUDES   = -I$(OPENCL_INCDIR) -I.
LIB        = -L$(OPENCL_LIBDIR) -lOpenCL
CXXFLAGS   = -fopenmp -O3


SOURCES_CPP =TestCalib.cpp
OBJECTS     =TestCalib.o  
EXECUTABLE  =calibrate

#g++ -I. -I/usr/local/cuda/include  -O3 tridag.cpp -lOpenCL


.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

default: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)


#PATH=.:$OCL_ROOT/bin:$PATH
#LD_LIBRARY_PATH=$OCL_ROOT/lib64:$OCL_ROOT/lib:$OCL_UTIL_ROOT/C/lib:$LD_LIBRARY_PATH
#
#CPLUS_INCLUDE_PATH=$OCL_ROOT/include:$OCL_UTIL_ROOT/C/common/inc:$OCL_UTIL_ROOT/OpenCL/common/inc:$OCL_UTIL_ROOT/shared/inc:$CPLUS_INCLUDE_PATH
