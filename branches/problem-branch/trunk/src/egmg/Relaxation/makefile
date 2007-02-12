Objects = LineRelaxation.o Jacobi.o GSLexicographic.o GSRedBlack.o ZebraLineGS.o LineGS.o LineJAC.o ZebraLine.o WHighJac.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
Jacobi.h GSLexicographic.h GSRedBlack.h LineRelaxation.h: Relaxation.h 
ZebraLineGS.h LineGS.h ZebraLine.h LineJAC.h: LineRelaxation.h
LineJAC.h: Jacobi.h
LineGS.h: GSLexicographic.h
ZebraLine.h ZebraLineGS.h: GSRedBlack.h
WHighJac.h: ../functions/WHighJacScalarProduct.h
Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
