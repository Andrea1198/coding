#
#=====================
# makefile_trial
#=====================
#

include ../make.inc

all: lib_hello

lib_hello:
	if test -d hello ; then (cd hello; $(MAKE) ) ; fi

clean wash:
	if test -d hello ; then (cd hello; $(MAKE) $@ ; fi
