Objects = TransferOperator.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

TransferOperator.h: ../Stencil/Stencil.h parameters.h

all: $(Objects)

clean:
	$(RM) $(Objects)
