Objects =
	
CycleType.h, FCycle.h, WCycle.h: Stencil.h
CycleType.h, FCycle.h, WCycle.h: ../general/parameters.h
VCycle.h: WCycle.h

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)