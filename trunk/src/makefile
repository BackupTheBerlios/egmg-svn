CXX = LANG=en g++
LD = LANG=en g++
export CXXFLAGS = -g -O3 -W -Wall -ansi -pedantic #-O3 -W -Wall -ansi -pedantic
LDFLAGS = $(CXXFLAGS)

Objects = main.o
Obj = ./functions/*.o ./Prolongation/*.o ./Relaxation/*.o ./Restriction/*.o ./Stencil/*.o ./Function/*.o
Exe = ./main.exe

all: TARGET = all
all: CycleType functions Prolongation Relaxation Restriction Stencil Function $(Objects)
	$(LD) $(LDFLAGS) $(Objects) $(Obj) -o $(Exe)
	
clean: TARGET = clean
clean: CycleType functions Prolongation Relaxation Restriction Stencil Function
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

.PHONY: CycleType
CycleType:
	$(MAKE) -C CycleType $(TARGET)
.PHONY: functions
functions:
	$(MAKE) -C functions $(TARGET)
.PHONY: Prolongation
Prolongation:
	$(MAKE) -C Prolongation $(TARGET)
.PHONY: Relaxation
Relaxation:
	$(MAKE) -C Relaxation $(TARGET)
.PHONY: Restriction
Restriction:
	$(MAKE) -C Restriction $(TARGET)
.PHONY: Stencil
Stencil:
	$(MAKE) -C Stencil $(TARGET)
.PHONY: Function
Function:
	$(MAKE) -C Function $(TARGET)

#always build main to keep dependency simple while testing
main.o: main.cpp CycleType functions Prolongation Relaxation Restriction Stencil Function
	$(CXX) $(CXXFLAGS) -c main.cpp
