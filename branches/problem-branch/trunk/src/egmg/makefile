all: TARGET = all
all: CycleType functions Prolongation Relaxation Restriction Stencil Function
	
clean: TARGET = clean
clean: CycleType functions Prolongation Relaxation Restriction Stencil Function

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
