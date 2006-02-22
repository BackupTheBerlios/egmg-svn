Objects = cycle.o error.o putbval.o putrhs.o residuum.o max_residuum.o twonorm_residuum.o output.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
cycle.h: ../Prolongation/Prolongation.h ../Relaxation/Relaxation.h ../Restriction/Restriction.h

cycle.h error.h f.h g.h zero.h putbval.h putrhs.h max_residuum.h residuum.h output.h twonorm_residuum.h: ../general/parameters.h
cycle.h max_residuum.h residuum.h twonorm_residuum.h: ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
