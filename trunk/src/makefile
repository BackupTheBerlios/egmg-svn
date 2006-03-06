CXX = LANG=en g++
LD = LANG=en g++
export CXXFLAGS = -O3 -W -Wall -ansi -pedantic
LDFLAGS = $(CXXFLAGS)

Objects = main.o
Obj = ./functions/*.o ./Prolongation/*.o ./Relaxation/*.o ./Restriction/*.o
Exe = ./main.exe

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

#always build main to keep dependency simple while testing
main.o: main.cpp CycleType functions Prolongation Relaxation Restriction Stencil
	$(CXX) $(CXXFLAGS) -c main.cpp


all: TARGET = all
all: CycleType functions Prolongation Relaxation Restriction Stencil $(Objects)
	$(LD) $(LDFLAGS) $(Objects) $(Obj) -o $(Exe)

clean: TARGET = clean
clean: CycleType functions Prolongation Relaxation Restriction Stencil
	$(RM) $(Exe) $(Objects) doxygen.log
	
debug: CXXFLAGS = -O0 -g -W -Wall -ansi -pedantic
debug: all
	
profile: CXXFLAGS = -O3 -pg -W -Wall -ansi -pedantic
profile: all
	$(Exe)
	gprof $(Exe) > gprof.out.txt
	$(RM) gmon.out

#do always a rebuild of the documentation if doc is called
doc:
	doxygen Multigrid.Doxyfile

run: all
	$(Exe)
