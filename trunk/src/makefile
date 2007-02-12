CXX = LANG=en g++
LD = LANG=en g++
export CXXFLAGS = -g -O3 -W -Wall -ansi -pedantic #-O3 -W -Wall -ansi -pedantic
LDFLAGS = $(CXXFLAGS)

Objects = main.o
Obj = ./egmg/functions/*.o ./egmg/Prolongation/*.o ./egmg/Relaxation/*.o ./egmg/Restriction/*.o ./egmg/Stencil/*.o ./egmg/Function/*.o
Exe = ./main.exe

all: TARGET = all
all: egmg $(Objects)
	$(LD) $(LDFLAGS) $(Objects) $(Obj) -o $(Exe)
	
clean: TARGET = clean
clean: egmg
	$(RM) $(Exe) $(Objects) doxygen.log

run: all
	$(Exe)
	
debug: CXXFLAGS = -O0 -g -W -Wall -ansi -pedantic -DDEBUG
debug: all
	
profile: CXXFLAGS = -O3 -pg -W -Wall -ansi -pedantic -fprofile-arcs -ftest-coverage
profile: all
	$(Exe)
#	gprof $(Exe) > gprof.out.txt
#	$(RM) gmon.out

#do always a rebuild of the documentation if doc is called
doc:
	doxygen Multigrid.Doxyfile

.PHONY: egmg
egmg:
	$(MAKE) -C egmg $(TARGET)

#always build main to keep dependency simple while testing
main.o: main.cpp egmg
	$(CXX) $(CXXFLAGS) -c main.cpp
