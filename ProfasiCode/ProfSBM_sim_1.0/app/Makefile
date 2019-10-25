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

# Makefile for ProFASi applications. 
# Last updated: 13 January 2012
#
VERSION ?= 1.4.9

#Source and lib directories for compilation. By default they 
#are set to use the "current" versions of the code in the 
#source directory, so that you can test the result of any 
#changes to the core package immediately. 
SrcDir?=../model
LIBDIR?=$(SrcDir)/lib
IncludeDir?=../include
BinDir?=../bin

CC?= g++
MPICC?=mpicxx
EX?=
MPIEX?=.mex
MAKE ?= make

IFLAGS := -I../model
ifeq ($(BUILDTYPE),debug)
EX:=.dbg
MPIEX:=.mdb
LIBprofasi := profasi_$(VERSION)_debug
LIBprofasi_utils := prf_utils_$(VERSION)_debug
LIBprofasi_par := profasi_$(VERSION)_par_debug
else
LIBprofasi := profasi_$(VERSION)
LIBprofasi_utils := prf_utils_$(VERSION)
LIBprofasi_par := profasi_$(VERSION)_par
endif

LIBS := -L$(LIBDIR) -l$(LIBprofasi) -l$(LIBprofasi_utils)
PLIBS := -l$(LIBprofasi_par)

ZLIB ?= -lz

include ../Makefile.compilers

EXEC := bin/pdb_slices$(EX) \
	bin/prf_energies$(EX) \
	bin/mimiqa$(EX) \
	bin/generate_contact_file$(EX) \
	bin/contact_order$(EX) \
	bin/prf_convert$(EX) \
	bin/extract_snapshot$(EX) \
	bin/extract_props$(EX) \
	bin/make_traj$(EX) \
	bin/minimize$(EX) \
	bin/regularize$(EX)
#	bin/Retrace$(EX)
SIM := bin/BasicMCRun$(EX) \
	bin/SimTempRun$(EX) \
	bin/SimAnnealRun$(EX)\
	bin/WLRun$(EX)
PARALLEL_EXEC := bin/ParTempRun$(MPIEX) \
				bin/HamRepExRun$(MPIEX) \
				$(SIM:%$(EX)=%$(MPIEX))
				
UTILS := bin/his1dmerge$(EX) bin/his2drender$(EX) \
	bin/prf_his1d$(EX) bin/prf_his2d$(EX)

MISC := bin/conf_to_gzip$(EX)
#bin/StAn$(EX)


serial := $(UTILS) $(EXEC) $(MISC) $(SIM)
all: $(serial)

parallel: $(PARALLEL_EXEC)

bin/%$(EX): %.cc 
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $< $(LIBS)

bin/%$(MPIEX) : %.cc 
	$(MPICC) -DPARALLEL $(IFLAGS) $(MPICCFLAGS) -o $@ $< $(LIBS) $(PLIBS)

bin/conf_to_gzip$(EX): bin/%$(EX) : ./%.cc 
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $< -L$(LIBDIR) -l$(LIBprofasi_utils) $(ZLIB)

#bin/StAn$(EX): bin/%$(EX) : ./%.cc ./%.hh
#	$(CC) $(CFLAGS) -o $@ $< ObsHandler.o$(EX) -L$(LIBDIR) -l$(LIBprofasi) -l$(LIBprofasi_utils)
clean:
	rm -f *.o$(EX) bin/*$(EX) bin/*$(MPIEX)
superclean:
	rm -f *.o *~ bin/*

OS := $(shell uname)
install:
ifeq ($(OS),Linux)
	for exefile in `ls bin/* `; do \
	install --backup -s -m 0755  ./"$${exefile}" $(BinDir)/ ; \
	done
else
	for exefile in `ls bin/* `; do \
	install -b -s -m 0755  ./"$${exefile}" $(BinDir)/ ; \
	done
endif
	
