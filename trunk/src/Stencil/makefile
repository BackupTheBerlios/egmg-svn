Objects =
	
Laplacian2D2.h Laplacian2D4.h Helmholtz2D2 Biharmonic2D2 ConDiff2D1 ConDiff2Dkap ConvectionDiffusion2D1: Stencil.h
Stencil.h: ../general/parameters.h

Galerkin.o: Galerkin.h Galerkin.cpp
	$(CXX) $(CXXFLAGS) -c Galerkin.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)
