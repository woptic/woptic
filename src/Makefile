### woptic/src/Makefile
###
###    Compile woptic executables and copy to ../bin
###
### Copyright 2013-2015 Elias Assmann
###

## include ‘make.sys’, if present
-include ../make.sys

####### Configuration variables #######

LDLIBS  ?= -llapack -lblas
W2W_DIR ?= $(WIENROOT)

#######################################

## Some defaults & derived variables
SHELL = /bin/bash
.SUFFIXES:

W2W_SRCDIR  := $(W2W_DIR)/SRC_w2w
W2W_LIBDIR  := $(W2W_SRCDIR)/lib
W2W_UTILOBJ := $(W2W_LIBDIR)/util_w2w.o
W2W_UTARGET := lib/util_w2w.o
UTILOBJ     := util_woptic.o


## Internal variables

# all executables (= targets)
EXE = $(EXE_W2WU) $(EXE_SIMPLE)

# executables (and objects) that use ‘w2w_util.F’
EXE_W2WU = compute_vr convert_vr kanalysis refine_tetra woptic_main	\
	   joinham inwopcheck

# executables (and objects) that use ‘woptic_util.f90’
EXE_UTIL = compute_vr convert_vr refine_tetra woptic_main joinham	\
	   inwopcheck

EXE_SIMPLE = obtain_dist

### Phony targets
all: $(EXE) logo.txt logo.ansi ../doc/logo.tex

clean:
	rm -f *.mod *.o *__genmod.{f90,mod}

distclean: clean
	rm -f $(EXE) logo.{ansi,txt}

.PHONY: all clean distclean

### Rules for executables
$(W2W_UTILOBJ): $(W2W_SRCDIR)/util_w2w.F $(W2W_SRCDIR)/Makefile
	make -C $(W2W_SRCDIR) $(W2W_UTARGET)

$(EXE_UTIL): %: %.o $(UTILOBJ)

$(EXE_W2WU): %: %.o $(W2W_UTILOBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
	mv -v $@ ../bin

$(EXE_SIMPLE): %: %.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
	mv -v $@ ../bin

W2WU_USERS := $(addsuffix .o,$(EXE_W2WU))
UTIL_USERS := $(addsuffix .o,$(EXE_UTIL))

refine_tetra: indexx.o

### Rules for objects
$(W2WU_USERS): $(W2W_UTILOBJ)
$(W2WU_USERS): override FFLAGS += -I$(W2W_LIBDIR)

$(UTIL_USERS): $(UTILOBJ)

%.o: %.[fF]90
	$(FC) $(FFLAGS) -c $<

### Rules for logo variants
logo.txt logo.ansi ../doc/logo.tex: logo.pl
	./logo.pl

## Time-stamp: <2015-11-09 14:24:50 assman@faepop36.tu-graz.ac.at>
