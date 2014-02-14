CXX        = g++
INCLUDES   = -I$(CUDA_INCDIR) -I. -I../../include
LIB        = -L$(OPENCL_LIBDIR) -lOpenCL
CXXFLAGS   = -O3 -fopenmp

SOURCES_CPP =VolCalibOrig.cpp
OBJECTS     =VolCalibOrig.o  
EXECUTABLE  =CalibVol

#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL


default: 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
