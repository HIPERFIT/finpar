include ../../setup.mk

SOURCES_CPP =TestCalib.cpp
OBJECTS     =TestCalib.o  
EXECUTABLE  =calibrate

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

default: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LIB)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
