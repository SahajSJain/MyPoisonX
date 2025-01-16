VERSION = 1.0
BOPT = O

EXE  = MyPoisonX
OBJDIR = ./obj

all: $(EXE)

CPP        = /lib/cpp
CPPFLAGS   = -D MPI

FC         = mpiifort
COMMFLAGS  = -O0 -fbounds-check -fcheck=all -g -fbacktrace  
#COMMFLAGS = -O3
CFLAGS     = ${COMMFLAGS} -J$(OBJDIR) -I$(OBJDIR)  # -J and -I options specify where to place and find module files
LFLAGS     =  

OBJS = $(OBJDIR)/MODULE.INPUTDATA.o     \
       $(OBJDIR)/MODULE.COEFFDATA.o     \
       $(OBJDIR)/MODULE.FIELDDATA.o     \
       $(OBJDIR)/MODULE.FLOWDATA.o      \
       $(OBJDIR)/MODULE.MPIDATA.o       \
       $(OBJDIR)/CODE.MAIN.o            \
       $(OBJDIR)/CODE.STARTUP.o         \
       $(OBJDIR)/CODE.MPISUBS.o         \
       $(OBJDIR)/CODE.SETUP_FIELD_VARIABLES.o \
       $(OBJDIR)/CODE.MEMORY.ALLOCATE.o \
       $(OBJDIR)/CODE.BOUNDARY_CONDITIONS.o \
       $(OBJDIR)/CODE.INITIALIZE.o      \
       $(OBJDIR)/CODE.WRITE_DUMP.o      \
       $(OBJDIR)/CALCULATOR.JACOBI.o   \
       $(OBJDIR)/CALCULATOR.MPITEST.o \
       $(OBJDIR)/CALCULATOR.GSSOR.o   \
       $(OBJDIR)/CALCULATOR.RESIDUAL.o 

# Create OBJDIR if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

%.c.f90: %.F90
	$(CPP) $(CPPFLAGS) $<  $@

$(OBJDIR)/%.o: %.c.f90 | $(OBJDIR)
	$(FC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJS)
	${FC} ${LFLAGS} -o ${EXE} $(OBJS)
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod MyPoisonX
