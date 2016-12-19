### woptic/Makefile
###
###    woptic main Makefile
###
### Copyright 2015 Elias Assmann
###

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
