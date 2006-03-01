Objects = LineRelaxation.o GSLexicographic.o GSRedBlack.o ZebraLineGS.o LineGS.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
GSLexicographic.h GSRedBlack.h LineRelaxation.h: Relaxation.h
LineGS.h: GSLexicographic.h
ZebraLineGS.h: GSRedBlack.h
LineGS.h ZebraLineGS.h: LineRelaxation.h
Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
