### woptic/Makefile
###
###    woptic main Makefile
###
### Copyright 2015 Elias Assmann
###
### $Id: Makefile 399 2015-06-03 20:28:10Z assmann $

SUBDIRS := src doc

.PHONY: all clean $(SUBDIRS) dist

all: src

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
		rm -f $$dir/Makefile.orig; \
	done

distclean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir distclean; \
	done

	rm -rf bin/compute_vr bin/convert_vr bin/inwopcheck \
	       bin/joinham bin/kanalysis bin/obtain_dist \
	       bin/refine_tetra bin/woptic_main


## Time-stamp: <2015-06-01 22:48:46 elias@hupuntu>
