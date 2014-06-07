include ../../setup.mk

ifndef ($(HAVE_GPU))
include ../../platform.mk
endif

INCLUDES  += -I../includeC
GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) -D GPU_DEV_ID=$(GPU_DEVICE_ID) \
             -D GPU_LOC_MEM=$(GPU_LOCAL_MEM) -D GPU_CONST_MEM=$(GPU_CONST_MEM) \
             -D GPU_REG_MEM=$(GPU_REG_MEM) -D GPU_GLB_MEM=$(GPU_DEVICE_MEM) \
             -D GPU_TILE=$(GPU_LOCAL_MEM_PER_TH) -D GPU_CORES=$(GPU_NUM_CORES)

SOURCES_CPP =GenPricing.cpp
HELPERS     =StructGPU.h ../includeC/Constants.h ../includeC/Optimizations.h ../includeC/ParseInput.h ../includeC/TimeHelper.h
OBJECTS     =GenPricing.o  
EXECUTABLE  =GenPricing

#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL
default: gpu

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

gpu: $(OBJECTS) $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

# The CPU version is in: ../Orig_COpenMP. This folder is GPU only!
cpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

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

