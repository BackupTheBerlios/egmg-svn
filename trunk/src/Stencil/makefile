Objects = Stencil.o Biharmonic2D2.o ConDiff2D1.o ConDiff2Dkap.o ConvectionDiffusion2D1.o Helmholtz2D2.o Laplacian2D2.o Laplacian2D4.o MSV2D4.o PreCalculatedStencil.o Galerkin.o StencilBoard.o 
	
Laplacian2D2.h Laplacian2D4.h Helmholtz2D2.h Biharmonic2D2.h ConDiff2D1.h ConDiff2Dkap.h ConvectionDiffusion2D1.h MSV2D4.h PreCalculatedStencil.h: Stencil.h ../Prolongation/Prolongation.h ../Restriction/Restriction.h
StencilBoard.h Stencil.h: ../general/parameters.h
PreCalculatedStencil.h: StencilBoard.h
Galerkin.h: PreCalculatedStencil.h
Galerkin.cpp: ../functions/expansion.h ../functions/generatePositionArrays.h

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)
