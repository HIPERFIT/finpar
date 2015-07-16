# The MAKE_DIR allows the preprocessor to include the path of the make file in the source code sent to the compiler.
# This is handled by a macro which is defined as an argument in the GPU_OPTS. 
# The C/C++ code can then use the variable called "MAKEFILE_PATH" without declaring/defining it.
MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

include ../../setup.mk

ifndef ($(HAVE_GPU))
  ifeq ($(OS),Darwin)
    include ../../platform.macbookpro.mk
  else
    include ../../platform.mk
  endif
endif

INCLUDES  += -I../includeC
GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) -D GPU_DEV_ID=$(GPU_DEVICE_ID) \
             -D GPU_LOC_MEM=$(GPU_LOCAL_MEM) -D GPU_CONST_MEM=$(GPU_CONST_MEM) \
             -D GPU_REG_MEM=$(GPU_REG_MEM) -D GPU_GLB_MEM=$(GPU_DEVICE_MEM) \
             -D GPU_TILE=$(GPU_LOCAL_MEM_PER_TH) -D GPU_CORES=$(GPU_NUM_CORES) \
			 -D CURR_DIR_PATH='"$(MAKE_DIR)"'

SOURCES_CPP =GenPricing.cpp
HELPERS     =StructGPU.h ../includeC/Constants.h ../includeC/Optimizations.h ../includeC/ParseInput.h
OBJECTS     =GenPricing.o  
EXECUTABLE  =GenPricing


ifeq ($(HAVE_GPU),1)
#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL
default: gpu
else
# cannot build this one... but proceed (return 0)
no-gpu:	
	@echo "$(HAVE_GPU) \n*** No GPU configured, cannot build `pwd` ***\n"
endif

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

gpu: $(OBJECTS) $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

# The CPU version is in: ../Orig_COpenMP. This folder is GPU only!
# cpu: $(OBJECTS)
# 	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

run_small: 
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_medium: 
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_large: 
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE) 2> Debug.txt


#default: $(OBJECTS)
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS) ./*.ptx ./Debug.txt
	@# clean nVidia compiler cache
	rm -rf $(HOME)/.nv/ComputeCache/*

