Objects = LineRelaxation.o Jacobi.o GSLexicographic.o GSRedBlack.o ZebraLineGS.o LineGS.o LineJAC.o ZebraLine.o WHighJac.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
Jacobi.h GSLexicographic.h GSRedBlack.h LineRelaxation.h SPAI.h: Relaxation.h 
LineGS.h: GSLexicographic.h
ZebraLineGS.h ZebraLine.h: GSRedBlack.h
LineJAC.h LineGS.h ZebraLineGS.h ZebraLine.h: LineRelaxation.h
LineJAC.h: Jacobi.h
WHighJac.h: ../functions/WHighJacScalarProduct.h
LineJAC.cpp ZebraLine.h: ../functions/residuum.h
Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
