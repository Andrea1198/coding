include ../../make.inc

DIR=baselib
LIBOBJS_LOC=$(LIBOBJS)

OBJS=

exec= libdyson_base.a

all : $(exec)

libdyson_base.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean : 
	- /bin/rm -f *.x *.o *.a *core *.F90 *.mod work.pc work.pcl

wash : clean

include make.depend
