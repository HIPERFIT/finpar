# If we are run from outside Hipermark itself, a number of environment
# variables will be missing.  Set these to hopefully-working values to
# support standalone compilation.  Users should not need to modify this.
HIPERMARK_BENCHMARK_LIB_DIR  ?= ../../lib/
HIPERMARK_IMPLEMENTATION_DIR ?= .

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
  OPENCL_LIBDIR     ?= $(OPENCL_ROOTDIR)/lib64
  OPENCL_INCDIR	    ?= $(OPENCL_ROOTDIR)/include
  CXX        = g++
  LIB        = -L$(OPENCL_LIBDIR) -lOpenCL
  CXXFLAGS   = -DENABLE_OPENMP -fopenmp -O3
  INCLUDES   = -I$(OPENCL_INCDIR) -I. -I../../include
endif

# Users should not modify this.
CXXFLAGS    += -DHIPERMARK_IMPLEMENTATION_DIR='"$(HIPERMARK_IMPLEMENTATION_DIR)"'
CXXFLAGS    += -DHIPERMARK_BENCHMARK_LIB_DIR='"$(HIPERMARK_BENCHMARK_LIB_DIR)"'
CXXFLAGS    += -DHIPERMARK_LIB_DIR='"$(HIPERMARK_LIB_DIR)"'
