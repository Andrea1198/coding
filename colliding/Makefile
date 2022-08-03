#===============================
# AGWX package
#===============================
#
MAKE=make

#
# manual
#
default:
	@echo
	@echo  "  *** colliding simulation code *** "
	@echo
	@echo  "  to configure, type:  ./configure [<options>]"
	@echo  "  to compile,   type:  make <target>"
	@echo
	@echo  "  Possible <target>'s are: "
	@echo
	@echo  "     all                all the executables "
	@echo  "     deps               update fortran90 dependencies "
	@echo  "     clean              remove executables and objects"
	@echo  "     clean_test         clean up the test examples"
	@echo  "     distclean          revert distribution to the original status"
	@echo

#
# MAIN target
#
all: libsrc

deps:
	if test -x ./config/makedeps.sh ; then ./config/makedeps.sh ; fi

#
# LIBS and MODULES
#
libextlibs:
	if test -d extlibs ; then \
        ( cd extlibs ; $(MAKE) ) ; fi

libsrc: cps

baselib_: libextlibs ctools_
	if test -d src/baselib ; then \
        ( cd src/baselib ; $(MAKE) ) ; fi

ctools_: libextlibs
	if test -d src/ctools ; then \
        ( cd src/ctools ; $(MAKE) ) ; fi

cps:  libextlibs baselib_ ctools_
	if test -d src/cps ; then \
        ( cd src/cps ; $(MAKE) ) ; fi
#
# CLEAN UP
#
clean:
	if test -d extlibs ;   then ( cd extlibs;   $(MAKE) clean ) ; fi
	if test -d src ;       then ( cd src;       $(MAKE) clean ) ; fi
	- /bin/rm -rf bin/*.x ./bin/iotk
	- /bin/rm -rf *.dat

clean_all: wash
distclean: wash
wash:  clean wash_extlibs
	- /bin/rm -rf make.inc ./config/configure.msg \
                ./config/config.log ./config/config.status \
                ./config/*.lineno \
                ./include/configure.h ./include/fft_defs.h \
                ./include/c_defs.h    ./include/f_defs.h\
                */dum1 */dum2
	- touch make.inc

wash_extlibs:
	if test -d extlibs ;   then ( cd extlibs;   $(MAKE) wash ) ; fi

