/*******************************************************************************
    PROFASI: Protein Folding and Aggregation Simulator, Version 1.5
    Copyright (C) (2012)  Anders Irback and Sandipan Mohanty
    Email: profasi@thep.lu.se
    Home Page: http://cbbp.thep.lu.se/activities/profasi/
    Version control (git) : https://trac.version.fz-juelich.de/PROFASI

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License
    (see PROFASI/gpl.txt).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
********************************************************************************/

#include "TrajSnapshot.hh"
#include "BinStream.hh"

using namespace prf_traj;
TrajSnapshot::TrajSnapshot() : tindex(0), ene(0), icyc(0), ncalls(0), seed(0)
{
}

TrajSnapshot::TrajSnapshot(const TrajSnapshot &t) : tindex(t.tindex),
ene(t.ene), icyc(t.icyc), ncalls(t.ncalls), seed(t.seed), popconf(t.popconf)
{
}

TrajSnapshot::~TrajSnapshot() {}

TrajSnapshot & TrajSnapshot::operator =(const TrajSnapshot &t)
                                       {
    if (this!=&t) {
        tindex=t.tindex;
        ene=t.ene;
        icyc=t.icyc;
        ncalls=t.ncalls;
        seed=t.seed;
        popconf=t.popconf;
    }
    return *this;
}

void TrajSnapshot::fill(BinStream &bst)
{
    bst>>icyc>>tindex>>ene>>seed>>ncalls;
    for (size_t i=0;i<popconf.size();++i) {
        bst>>popconf[i];
    }
}
