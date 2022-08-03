include ../make.inc

all: lib_base lib_gf

lib_ctools : 
	if test -d ctools ; 	then ( cd ctools ; 	$(MAKE) ) ; fi

lib_base: ctools
	if test -d baselib ; 	then ( cd baselib ; $(MAKE) ) ; fi

lib_gf: lib_base
	if test -d gf ; 		then ( cd gf ; 		$(MAKE) ) ; fi

clean wash:
	if test -d ctools ; 	then ( cd ctools;   $(MAKE) $@ ) ; fi
	if test -d baselib ; 	then ( cd baselib; 	$(MAKE) $@ ) ; fi
	if test -d gf ; 		then ( cd gf;		$(MAKE) $@ ) ; fi
