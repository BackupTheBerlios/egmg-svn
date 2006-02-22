Objects = Injection.o Full_Weighting.o Half_Weighting.o seven_point_restr.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
Injection.h Full_Weighting.h Half_Weighting.h seven_point_restr.h: Restriction.h

Restriction.h: ../general/parameters.h ../general/Transfer_Operator.h

all: $(Objects)

clean:
	$(RM) $(Objects)
