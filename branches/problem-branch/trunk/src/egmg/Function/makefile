Objects = Function.o Zero.o Constant.o TestFunction1.o TestRightSide1.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

all: $(Objects)

clean:
	$(RM) $(Objects)
