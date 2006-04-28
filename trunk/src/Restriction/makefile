Objects = Injection.o FullWeighting.o HalfWeighting.o SevenPointRestriction.o TranspProlong.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
Injection.h FullWeighting.h HalfWeighting.h SevenPointRestriction.h TranspProlong.h: Restriction.h

Restriction.h: ../general/parameters.h ../general/TransferOperator.h ../Stencil/Stencil.h ../Prolongation/Prolongation.h

all: $(Objects)

clean:
	$(RM) $(Objects)
