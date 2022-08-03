include ../../make.inc


obj=\
	errore.o \
	qties.o \
	update.o

all: lib_collision.a cse.x

cse.x: lib_make.a
	$(F90COMP) -o $@ cse.o $(obj)
	( ln -f $@ ../../bin/cse.x )

lib_collision.a:
	gfortran -c *.f90

clear:
	rm $(obj)
	rm main.x
