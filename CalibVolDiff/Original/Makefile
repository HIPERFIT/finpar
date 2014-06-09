CXX        = g++
INCLUDES   = -I$(CUDA_INCDIR) -I. -I../../include
LIB        = -L$(OPENCL_LIBDIR) -lOpenCL
CXXFLAGS   = -O3 -ffast-math  -fopenmp

SOURCES_CPP =VolCalibOrig.cpp
OBJECTS     =VolCalibOrig.o  
EXECUTABLE  =CalibVol

#g++ -I. -I/usr/local/cuda/include  -O3 NordeaVect.cpp -lOpenCL


default: 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)

run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_medium: $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE) 2> Debug.txt

run_large: $(EXECUTABLE)
	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE) 2> Debug.txt

clean:
	rm -f Debug.txt $(EXECUTABLE) $(OBJECTS)
