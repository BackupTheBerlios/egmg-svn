Objects =
	
Laplacian2D2.h Laplacian2D4.h Helmholtz2D2.h Biharmonic2D2.h ConDiff2D1.h ConDiff2Dkap.h ConvectionDiffusion2D1.h MSV2D4.h: Stencil.h ../Prolongation/Prolongation.h ../Restriction/Restriction.h
Stencil.h: ../general/parameters.h

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)
