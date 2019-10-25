#
#   PROFASI: Protein Folding and Aggregation Simulator, Version 1.5
#   Copyright (C) (2012)  Anders Irback and Sandipan Mohanty
#   Email: profasi@thep.lu.se
#   Home Page: http://cbbp.thep.lu.se/activities/profasi/
#   Version control (git) : https://trac.version.fz-juelich.de/PROFASI

#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License
#   (see PROFASI/gpl.txt).

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#

# Makefile for profasi/model.  
# Last revision: 13 January 2012

VERSION ?= 1.4.9

BUILDTYPE ?= optimized
LibDir ?= ../lib
IncludeDir ?= ../include
SrcDir?=.

MAKE ?= make
CC ?= g++
MPICC ?= mpicxx
MPIOPTS ?= 

AR ?= ar
AROPTS ?= -r 
RANLIB ?= ranlib

include ../Makefile.compilers

ORD ?= obj

ifeq ($(BUILDTYPE),debug)
EX:=.dbg
LIBprofasi := lib/libprofasi_$(VERSION)_debug.a
LIBprofasi_utils := lib/libprf_utils_$(VERSION)_debug.a
LIBprofasi_par := lib/libprofasi_$(VERSION)_par_debug.a
else
LIBprofasi := lib/libprofasi_$(VERSION).a
LIBprofasi_utils := lib/libprf_utils_$(VERSION).a
LIBprofasi_par := lib/libprofasi_$(VERSION)_par.a
endif

MODULES := Aux Elements AA EndGroups Energy Energy/Extras Updates Algorithms Observables
include $(patsubst %,%/module.mk,$(MODULES))
OBJ := $(patsubst %.cc,$(ORD)/%.o,$(SRC))
UOBJ := $(patsubst %.cc,$(ORD)/%.o,$(USRC))
POBJ := $(patsubst %.cc,$(ORD)/%.o,$(PSRC))

all: $(LIBprofasi_utils) $(LIBprofasi) 
parallel: $(LIBprofasi_par)

$(LIBprofasi_utils): $(UOBJ)
	$(AR) $(AROPTS) $@ $^
	$(RANLIB) $@
$(LIBprofasi): $(OBJ)
	$(AR) $(AROPTS) $@ $^
	$(RANLIB) $@
$(LIBprofasi_par): $(POBJ)
	$(AR) $(AROPTS) $@ $^
	$(RANLIB) $@

$(ORD)/%.o : %.cc
	$(CC) -c $(CFLAGS) $< -o $@
$(POBJ): $(ORD)/%.o : %.cc
	$(MPICC) -c $(MPICCFLAGS) $< -o $@
clean: 
	rm -f $(LIBprofasi) $(LIBprofasi_utils) $(LIBprofasi_par)
	for i in $(MODULES); do \
		rm -f $(ORD)/$${i}/*.o ; \
	done
	
superclean: 
	rm -f $(LIBprofasi) $(LIBprofasi_utils) $(LIBprofasi_par)
	for i in $(MODULES); do \
		rm -f $${i}/*~ $${i}/*.backup $${i}/*.d; \
		rm -f $(ORD)/$${i}/*.o ; \
	done

OS := $(shell uname)
install:
ifeq ($(OS),Linux)
	install --backup --mode=0644 $(LIBprofasi_utils) $(LibDir)/
	install --backup --mode=0644 $(LIBprofasi) $(LibDir)/
	install --backup --mode=0644 $(LIBprofasi_par) $(LibDir)/
	for i in $(MODULES); do \
	mkdir -p $(IncludeDir)/$${i} ;\
	install --backup --mode=0644  $(SrcDir)/$${i}/*.hh $(IncludeDir)/$${i} ;\
	done
	install --backup --mode=0644  $(SrcDir)/Aux/*.tcc  $(IncludeDir)/Aux ;\
	install --backup --mode=0644  $(SrcDir)/Elements/*.i.cc $(IncludeDir)/Elements
	install --backup --mode=0644  $(SrcDir)/*.hh $(IncludeDir)
else
	install -b -m 0644 $(LIBprofasi_utils) $(LibDir)/
	install -b -m 0644 $(LIBprofasi) $(LibDir)/
	install -b -m 0644 $(LIBprofasi_par) $(LibDir)/
	for i in $(MODULES); do \
	mkdir -p $(IncludeDir)/$${i} ;\
	install -b -m 0644  $(SrcDir)/$${i}/*.hh $(IncludeDir)/$${i} ;\
	done
	install -b -m 0644  $(SrcDir)/Aux/*.tcc  $(IncludeDir)/Aux ;\
	install -b -m 0644  $(SrcDir)/Elements/*.i.cc $(IncludeDir)/Elements
	install -b -m 0644  $(SrcDir)/*.hh $(IncludeDir)
endif
include $(OBJ:.o=.d) $(UOBJ:.o=.d)
$(ORD)/%.d: %.cc 
	../bin/depend.sh $(ORD)/`dirname $*.cc` $*.cc > $@

