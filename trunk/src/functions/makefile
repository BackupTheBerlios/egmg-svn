Objects = cycle.o error.o putbval.o putrhs.o residuum.o max_residuum.o twonorm_residuum.o output.o

cycle.o: cycle.h cycle.cpp
	$(CXX) $(CXXFLAGS) -c cycle.cpp
	
error.o: error.h error.cpp
	$(CXX) $(CXXFLAGS) -c error.cpp
	
putbval.o: putbval.h putbval.cpp
	$(CXX) $(CXXFLAGS) -c putbval.cpp

putrhs.o: putrhs.h putrhs.cpp
	$(CXX) $(CXXFLAGS) -c putrhs.cpp
	
residuum.o: residuum.h residuum.cpp
	$(CXX) $(CXXFLAGS) -c residuum.cpp

max_residuum.o: max_residuum.h max_residuum.cpp
	$(CXX) $(CXXFLAGS) -c max_residuum.cpp

twonorm_residuum.o: twonorm_residuum.cpp  twonorm_residuum.h
	$(CXX) $(CXXFLAGS) -c twonorm_residuum.cpp

output.o: output.cpp output.h
	$(CXX) $(CXXFLAGS) -c output.cpp

cycle.h: ../Prolongation/Prolongation.h ../Relaxation/Relaxation.h ../Restriction/Restriction.h

cycle.h error.h f.h g.h zero.h putbval.h putrhs.h max_residuum.h residuum.h output.h twonorm_residuum.h: ../general/parameters.h
cycle.h max_residuum.h residuum.h twonorm_residuum.h: ../Stencil/Stencil.h

all: $(Objects)

clean:
	$(RM) $(Objects)
