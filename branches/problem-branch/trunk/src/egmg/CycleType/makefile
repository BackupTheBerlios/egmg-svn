Objects =

CycleType.h: Stencil.h ../general/parameters.h
FCycle.h, WCycle.h: CycleType.h
VCycle.h: WCycle.h

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)