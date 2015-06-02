include ../../setup.mk

ifndef ($(HAVE_GPU))
include ../../platform.mk
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
gpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

cpu: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)


run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_medium: $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_large: $(EXECUTABLE)
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE) 2> Debug.txt

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

