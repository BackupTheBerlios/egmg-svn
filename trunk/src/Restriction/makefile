Objects = Injection.o Full_Weighting.o Half_Weighting.o seven_point_restr.o

Injection.o: Injection.h Injection.cpp
	$(CXX) $(CXXFLAGS) -c Injection.cpp
Full_Weighting.o: Full_Weighting.h Full_Weighting.cpp
	$(CXX) $(CXXFLAGS) -c Full_Weighting.cpp
Half_Weighting.o: Half_Weighting.h Half_Weighting.cpp
	$(CXX) $(CXXFLAGS) -c Half_Weighting.cpp
seven_point_restr.o: seven_point_restr.h seven_point_restr.cpp
	$(CXX) $(CXXFLAGS) -c seven_point_restr.cpp
	
Injection.h Full_Weighting.h Half_Weighting.h seven_point_restr.h: Restriction.h

Restriction.h: ../general/parameters.h ../general/Transfer_Operator.h

all: $(Objects)

clean:
	$(RM) $(Objects)
