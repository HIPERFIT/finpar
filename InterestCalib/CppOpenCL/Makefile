MAKE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

include ../../setup.mk

ifndef ($(HAVE_GPU))
include ../../platform.mk
endif

INCLUDES  += -I../includeC
GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) -D GPU_DEV_ID=$(GPU_DEVICE_ID) \
             -D GPU_LOC_MEM=$(GPU_LOCAL_MEM) -D GPU_CONST_MEM=$(GPU_CONST_MEM) \
             -D GPU_REG_MEM=$(GPU_REG_MEM) -D GPU_GLB_MEM=$(GPU_DEVICE_MEM) \
             -D GPU_TILE=$(GPU_LOCAL_MEM_PER_TH) -D GPU_CORES=$(GPU_NUM_CORES) \
			 -D CURR_DIR_PATH='"$(MAKE_DIR)"'

SOURCES_CPP =SwapCalib.cpp
HELPERS     =../includeC/ParseInput.h
OBJECTS     =SwapCalib.o  
EXECUTABLE  =SwapCalib


default: gpu

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

# The GPU version is in: ../CppAndGPU. This folder is CPU only!
gpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

cpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)


run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_medium: $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_large: $(EXECUTABLE)
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE) 2> Debug.txt

clean:
	rm -f $(EXECUTABLE) $(OBJECTS) ./*.ptx SrcCL/*.ptx ./Debug.txt
	@# clean nVidia compiler cache
	rm -rf $(HOME)/.nv/ComputeCache/*

