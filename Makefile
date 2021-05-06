NAME=test

# use module_timer.f90
TIMER=on

#module switch PrgEnv-cray PrgEnv-intel
#

FC=ifort
#FC=ftn
#FC=gfortran
FFLAGS=-O3  -axMIC-AVX512 -mcmodel=large -qopenmp #-parallel
#FFLAGS=-O3 -axMIC-AVX512 -qopenmp -mcmodel=medium -shared-intel -fpic -dynamic
#FFLAGS=-O0   -CB -fpe0 -traceback -g -check uninit #-check noarg_temp_created #-heap-arrays -warn
FFLAGS2=$(FFLAGS)
MPIC=mpirun -np

ifeq ("$(TIMER)","on")
  FFLAGS+=-D_TIMER_ -D_TIMER_NOMPI_
endif
FFLAGS2=$(FFLAGS)
AR="ar -r"
MOD_DIR = ./
LIB_DIR = ./
FFLAGS+=-cpp


OBJS = parameter.o main.o set_all.o poisson.o helmholtz.o \
       check.o ibm.o file_make.o rg.o

ifeq ("$(TIMER)","on")
OBJS +=module_timer.o
endif

#.SUFFIXES : .f90
%.o : %.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c  $*.f90
#.SUFFIXES : .f
%.o : %.f
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c  $*.f
#.SUFFIXES : .F
%.o : %.F
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c  $*.F

$(NAME).x :  $(OBJS)  
	$(FC) $(FFLAGS2) -o  $(NAME).x $(OBJS) -mkl #-lfftw3 -lfftw3_threads

clean :
	rm -f $(OBJS)  *~ core*  *.mod *.o fort.* *.lst *.L *.a $(NAME).x

run  : $(NAME).x
	./$(NAME).x

ifeq ("$(TIMER)","on")
main.o:module_timer.o
march.o:module_timer.o
endif
