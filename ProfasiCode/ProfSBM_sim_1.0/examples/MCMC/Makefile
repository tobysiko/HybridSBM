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
#

CC= g++

VERSION = 1.1.6
EXECSUFFIX=ex

ProFASiDir=$(HOME)/profasi
LIBDIR=$(ProFASiDir)/lib
INCLUDEDIR=$(ProFASiDir)/include

#CFLAGS= $(NEED_MPI_CPP_BINDING_OPT) -g3 -Wall -I$(SrcDir)  
CFLAGS= -O3 -Wall -I$(INCLUDEDIR)  
RM = rm
Make = make

LIBprofasi=-lprofasi_$(VERSION)
LIBprofasi_par=-lprofasi_$(VERSION)_par
EXEC =  
RENDER = mcmc.$(EXECSUFFIX)

$(RENDER) : %.ex : %.cc 
	$(CC) $(CFLAGS) $< -o $@ -L$(LIBDIR) $(LIBprofasi)


clean: 
	rm -f *.o *~ *.$(EXECSUFFIX) *.$(MPIEXECSUFFIX)



