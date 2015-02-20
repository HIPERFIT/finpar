OS=$(shell uname -s)

ifeq ($(OS),Darwin)
  OPENCL_ROOTDIR    ?= /System/Library/Frameworks/OpenCL.framework
  OPENCL_LIBDIR     := $(OPENCL_ROOTDIR)/Libraries
  OPENCL_INCDIR	    ?= $(OPENCL_ROOTDIR)/Headers
  CXX        = g++-4.9
  CXXFLAGS   = 
  LIB        = -framework OpenCL -L$(OPENCL_LIBDIR) -Wl,-x -m64
  CXXFLAGS   = -fopenmp -O3 -fno-rtti
else
  OPENCL_ROOTDIR    ?= /usr/local/cuda
  OPENCL_LIBDIR     ?= $(OPENCL_ROOTDIR)/lib64
  OPENCL_INCDIR	    ?= $(OPENCL_ROOTDIR)/include
  CXX        = g++
  LIB        = -L$(OPENCL_LIBDIR) -lOpenCL
  CXXFLAGS   = -fopenmp -O3
endif

INCLUDES   = -I$(OPENCL_INCDIR) -I. -I../../include
