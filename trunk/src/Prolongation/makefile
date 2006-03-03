Objects = BicubicInterpolation.o BilinearInterpolation.o SevenPointInterpolation.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
BilinearInterpolation.h: Prolongation.h
SevenPointInterpolation.h: Prolongation.h
BicubicInterpolation.h: Prolongation.h

Prolongation.h: ../general/parameters.h

all: $(Objects)

clean:
	$(RM) $(Objects)
