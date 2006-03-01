Objects = GSLexicographic.o GSRedBlack.o ZebraLineGS.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
GSLexicographic.h GSRedBlack.h ZebraLineGS.h: Relaxation.h
ZebraLineGS.h: GSRedBlack.h

Relaxation.h: ../general/parameters.h ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
