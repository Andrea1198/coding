# install/make_iotk.inc.  Generated from make_iotk.inc.in by configure.

# -----------------------------
# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

# -----------------------------
# compiler definitions


CC                  = gcc
CFLAGS              = -O3 $(DFLAGS) $(IFLAGS)
CPP                 = cpp
CPPFLAGS            = -P -traditional $(DFLAGS) $(IFLAGS)
F90                 = gfortran
MPIF90              = gfortran
F90FLAGS            = $(FFLAGS) -x f95-cpp-input $(FDFLAGS) $(IFLAGS) 
F77                 = gfortran
FFLAGS              = -O3 -g
FFLAGS_NOOPT        = -O0 -g
LD                  = gfortran
LDFLAGS             = -g -pthread
AR                  = ar
ARFLAGS             = ruv
RANLIB              = ranlib


#------------------------------
# specific definitions

DFLAGS         =  -D__FFTW
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)
IFLAGS         = $(AUTO_INCLUDE)
LIBOBJS        = $(AUTO_LIBS)


# -----------------------------
# libraries

TOPDIR              = /e/coding/makefile_trial
#
LIB_MPI             = 
LIB_MASS            = 
#
LIBS                = $(LIB_MPI) $(LIB_MASS)

