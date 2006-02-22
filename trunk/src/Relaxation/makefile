Objects = GS_Lex.o GS_RB.o lineGS.o lineJAC.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
GS_Lex.h GS_RB.h lineGS.h lineJAC.h: Relaxation.h

Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
