Objects = GS_Lex.o GS_RB.o lineGS.o lineJAC.o

GS_Lex.o: GS_Lex.h GS_Lex.cpp
	$(CXX) $(CXXFLAGS) -c GS_Lex.cpp
GS_RB.o: GS_RB.h GS_RB.cpp
	$(CXX) $(CXXFLAGS) -c GS_RB.cpp
lineGS.o: lineGS.h lineGS.cpp
	$(CXX) $(CXXFLAGS) -c lineGS.cpp
lineJAC.o: lineJAC.h lineJAC.cpp
	$(CXX) $(CXXFLAGS) -c lineJAC.cpp
	
GS_Lex.h GS_RB.h lineGS.h lineJAC.h: Relaxation.h

Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
