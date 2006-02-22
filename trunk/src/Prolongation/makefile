Objects = Bilinear_interpolation.o SevenPoint_interpolation.o Bicubic_interpolation.o Dendy_interpolation.o DeZeeuw_interpolation.o

%.o: %.h %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp
	
Bilinear_interpolation.h: Prolongation.h
SevenPoint_interpolation.h: Prolongation.h
Bicubic_interpolation.h: Prolongation.h
Dendy_interpolation.h: Prolongation.h
DeZeeuw_interpolation.h: Prolongation.h

Prolongation.h: ../general/parameters.h

all: $(Objects)

clean:
	$(RM) $(Objects)
