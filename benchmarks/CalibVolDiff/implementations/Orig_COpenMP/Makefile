include ../../setup.mk

ifndef ($(HAVE_GPU))
include ../../platform.mk
endif

INCLUDES    += -I../includeC
GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) 

SOURCES_CPP =VolCalibOrig.cpp
HELPERS     =../includeC/ParseInput.h
OBJECTS     =VolCalibOrig.o  
EXECUTABLE  =VolCalib


default: cpu

.cpp.o: $(SOURCES_CPP) $(HELPERS)
	$(CXX) $(CXXFLAGS) $(GPU_OPTS) $(INCLUDES) -c -o $@ $<

# The GPU version is in: ../VectAll and ../VectOuters. This folder is CPU only!

cpu: $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS)

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)

