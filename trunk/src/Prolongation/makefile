Objects = Bilinear_interpolation.o SevenPoint_interpolation.o Bicubic_interpolation.o Dendy_interpolation.o DeZeeuw_interpolation.o

Bilinear_interpolation.o: Bilinear_interpolation.h Bilinear_interpolation.cpp
	$(CXX) $(CXXFLAGS) -c Bilinear_interpolation.cpp
	
SevenPoint_interpolation.o: SevenPoint_interpolation.h SevenPoint_interpolation.cpp
	$(CXX) $(CXXFLAGS) -c SevenPoint_interpolation.cpp
	
Bicubic_interpolation.o: Bicubic_interpolation.h Bicubic_interpolation.cpp
	$(CXX) $(CXXFLAGS) -c Bicubic_interpolation.cpp
	
Dendy_interpolation.o: Dendy_interpolation.h Dendy_interpolation.cpp
	$(CXX) $(CXXFLAGS) -c Dendy_interpolation.cpp
	
DeZeeuw_interpolation.o: DeZeeuw_interpolation.h DeZeeuw_interpolation.cpp
	$(CXX) $(CXXFLAGS) -c DeZeeuw_interpolation.cpp
	
Bilinear_interpolation.h: Prolongation.h
SevenPoint_interpolation.h: Prolongation.h
Bicubic_interpolation.h: Prolongation.h
Dendy_interpolation.h: Prolongation.h
DeZeeuw_interpolation.h: Prolongation.h

Prolongation.h: ../general/parameters.h

all: $(Objects)

clean:
	$(RM) $(Objects)
