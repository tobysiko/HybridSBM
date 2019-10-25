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

# Makefile for programs using the profasi library.

CC= g++

VERSION = 1.4.0
EXECSUFFIX=ex

ProFASiDir=../..
LIBDIR=$(ProFASiDir)/lib
INCLUDEDIR=$(ProFASiDir)/include

# Recommended compiler options for Intel icc compiler
ifeq ($(CC),icc)
CFLAGS=-O3 -fno-exceptions -w1 -I$(INCLUDEDIR)
MPICCFLAGS= -O3 -w1 -I$(INCLUDEDIR)
endif

# Recommended compiler options for GNU g++ compiler
ifeq ($(CC),g++)
CFLAGS=-g3  -Wall -I$(INCLUDEDIR)
MPICCFLAGS= $(CFLAGS)
endif

# Recommended compiler options for Portland Groups compiler suite
ifeq ($(CC),pgCC)
CFLAGS=-fast --no_exceptions -Minline=levels:10  -I$(INCLUDEDIR)
MPICCFLAGS=-fast -Minline=levels:10 -I(INCLUDEDIR)
endif

# Recommended compiler options for IBM xlC compiler
ifeq ($(CC),xlC)
CFLAGS=-O3 -qhot -I$(INCLUDEDIR)
MPICCFLAGS=$(CFLAGS)
endif

RM = rm
Make = make

LIBprofasi=-lprf_utils_$(VERSION)
RENDER = his_example.$(EXECSUFFIX) 
all: $(RENDER)

$(RENDER) : %.ex : %.cc 
	$(CC) $(CFLAGS) $< -o $@ -L$(LIBDIR) $(LIBprofasi)


clean: 
	rm -f *.o *~ *.$(EXECSUFFIX)



