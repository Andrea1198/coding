MAKE=make

default:
		@echo 
		@echo " *** Makefile trial code *** "
		@echo
		@echo " to configure type ./configure "
		@echo " to compile type make <target> "
		@echo
		@echo " <target>'s :"
		@echo "		all				all executables"
		@echo "		deps			update fortran90 dependencies"
		@echo "		clean			remove executables and objects"
		@echo "		clean_test		clean up test examples"
		@echo "		distclean		revert distribution to the original status"

all : 
	cd src ; $(MAKE) 

deps:
	if test -x ./config/makedeps.sh; then ./config/makedeps.sh ; fi

clean:
	if test -d src ; then ( cd src; $(MAKE) clean ) ; fi

clean_all: wash
distclean: wash


wash: clean
