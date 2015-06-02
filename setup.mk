OS=$(shell uname -s)

ifeq ($(OS),Darwin)
  # OpenMP is not supported by clang and gcc-4.9 does not work
  # well with OpenCL on Mac, thus, we don't define ENABLE_OPENMP
  CXX        = clang++
  CXXFLAGS   = -Wall -W -O3
  LIB        = -framework OpenCL
  INCLUDES   = -I. -I../../include
else
  ENABLE_OPENMP     = 1
  OPENCL_ROOTDIR    ?= /usr/local/cuda
  OPENCL_LIBDIR     := $(OPENCL_ROOTDIR)/lib64
  OPENCL_INCDIR	    ?= $(OPENCL_ROOTDIR)/include
  CXX        = g++
  LIB        = -L$(OPENCL_LIBDIR) -lOpenCL
  CXXFLAGS   = -DENABLE_OPENMP -fopenmp -O3
  INCLUDES   = -I$(OPENCL_INCDIR) -I. -I../../include
endif

