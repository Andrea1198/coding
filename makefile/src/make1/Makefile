include ../../make.inc

SOURCES := $(shell find $(SOURCEDIR) -name '*.f90')

obj=hello.o

all: lib_make.a main.x

main.x: lib_make.a
	$(F90COMP) -o $@ main.o $(obj)
	( ln -f main.x ../../bin/main1.x )

lib_make.a:
	$(F90COMP) -c *.f90

clear:
	rm $(obj)
	rm main.x
