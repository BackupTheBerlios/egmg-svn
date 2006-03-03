Objects = LineRelaxation.o Jacobi.o GSLexicographic.o GSRedBlack.o ZebraLineGS.o LineGS.o LineJAC.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
Jacobi.h GSLexicographic.h GSRedBlack.h LineRelaxation.h: Relaxation.h
LineGS.h: GSLexicographic.h
ZebraLineGS.h: GSRedBlack.h
LineJAC.h LineGS.h ZebraLineGS.h: LineRelaxation.h
LineJAC.h: Jacobi.h
LineJAC.cpp: ../functions/residuum.h
Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
