Objects = BicubicInterpolation.o BilinearInterpolation.o DendyInterpolation.o DeZeeuwInterpolation.o SevenPointInterpolation.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
BilinearInterpolation.h: Prolongation.h
SevenPointInterpolation.h: Prolongation.h
BicubicInterpolation.h: Prolongation.h
DendyInterpolation.h: Prolongation.h
DeZeeuwInterpolation.h: Prolongation.h

Prolongation.h: ../general/parameters.h

all: $(Objects)

clean:
	$(RM) $(Objects)
