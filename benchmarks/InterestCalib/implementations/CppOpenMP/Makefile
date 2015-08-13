include setup.mk

ifndef ($(HAVE_GPU))
include platform.mk
endif

INCLUDES    += -I../includeC
GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP)
#-fno-unsafe-math-optimizations -fno-finite-math-only 
#-fmath-errno -ftrapping-math -fno-rounding-math -fsignaling-nans

SOURCES_CPP =SwapCalib.cpp
HELPERS     =../includeC/ParseInput.h
OBJECTS     =SwapCalib.o  
EXECUTABLE  =SwapCalib

default: cpu

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

# The GPU version is in: ../CppAndGPU. This folder is CPU only!
cpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)
clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

