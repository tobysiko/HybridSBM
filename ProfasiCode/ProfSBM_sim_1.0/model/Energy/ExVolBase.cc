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

#include "ExVolBase.hh"

namespace prf
{
    ExVolBase::ExVolBase()
    {
        NA=5;
        ksa=0.1;
        cut=4.3;
        sigsa.resize(NA,0);
        sigsa[0]=1.00;//H
        sigsa[1]=1.75;//C
        sigsa[2]=1.53;//N
        sigsa[3]=1.42;//O
        sigsa[4]=1.77;//S
        cut2=cut*cut;
    }

    ExVolBase::~ExVolBase() {}
}
