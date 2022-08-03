f90comp = gfortran

mpicomp = mpif90

flags = 

AR = ar
ARFLAGS = ruv


objects = \
hello.o \
kinds.o \

all : main.x lib_hello.a

main.x : main.o lib_hello.a
	$(f90comp) $(flags) -o $@ main.o $(OBJECTS)
	
lib_hello.a : $(objects)
	$(f90comp) -c *.f90
	$(AR) $(ARFLAGS) $@ $?
	echo $@

clean:
	rm main.x $(objects)

include make.depend
