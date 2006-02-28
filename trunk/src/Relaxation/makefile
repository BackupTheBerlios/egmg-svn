Objects = GSLexicographic.o GSRedBlack.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
GSLexicographic.h GSRedBlack.h lineGS.h lineJAC.h: Relaxation.h

Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
