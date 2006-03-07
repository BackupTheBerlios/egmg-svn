Objects = cycle.o error.o putbval.o putrhs.o residuum.o maxResiduum.o twonormResiduum.o output.o directSolver.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
cycle.h: ../Prolongation/Prolongation.h ../Relaxation/Relaxation.h ../Restriction/Restriction.h
cycle.cpp: residuum.h directSolver.h
cycle.h error.h f.h g.h zero.h putbval.h putrhs.h maxResiduum.h residuum.h output.h twonormResiduum.h directSolver.h: ../general/parameters.h
cycle.h maxResiduum.h residuum.h twonormResiduum.h directSolver.h: ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
