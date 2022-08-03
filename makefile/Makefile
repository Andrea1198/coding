directories = make1 make2

MAKE = make

TOPDIR = `pwd`

all : make1 make2

make1 : 
	( cd $(TOPDIR)/src/$@ ; $(MAKE) )

make2 : 
	( cd $(TOPDIR)/src/$@ ; $(MAKE) )

clear : 
	for dir in $(directories) ; do \
		( cd $(TOPDIR)/src/$$dir ; $(MAKE) $@ ) \
	done
	rm ./bin/*

clear_all : clear
	
