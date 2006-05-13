Objects = Galerkin.o StencilBoard.o 
	
Laplacian2D2.h Laplacian2D4.h Helmholtz2D2.h Biharmonic2D2.h ConDiff2D1.h ConDiff2Dkap.h ConvectionDiffusion2D1.h MSV2D4.h Galerkin.h: Stencil.h ../Prolongation/Prolongation.h ../Restriction/Restriction.h
StencilBoard.h Stencil.h: ../general/parameters.h
Galerkin.h: StencilBoard.h ../Prolongation/Prolongation.h ../Restriction/Restriction.h
Galerkin.cpp: ../functions/expansion.h ../functions/generatePositionArrays.h

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)
