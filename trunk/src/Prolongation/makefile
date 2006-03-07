Objects = BicubicInterpolation.o BilinearInterpolation.o SevenPointInterpolation.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

BicubicInterpolation.h BilinearInterpolation.h SevenPointInterpolation.h: Prolongation.h

Prolongation.h: ../general/parameters.h ../general/TransferOperator.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
