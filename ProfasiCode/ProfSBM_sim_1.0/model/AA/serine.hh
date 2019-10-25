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

#ifndef Serine_HH
#define Serine_HH
#include "AminoAcid.hh"
#include "../Elements/ATetGroup.hh"
#include "../Elements/DihedralGroup.hh"

namespace prf
{
//! Serine
    /**
     * Specification of serine side chain. The OH group is represented by a
     * DihedralGroup, with just one out-going bond.
     * \ingroup aminoacids
     */

    class Serine:public AminoAcid
    {

    public:
        Serine();
        ~Serine();
        void Initialize();

    private:
        ATetGroup g0;
        DihedralGroup g1;
    };
}

#endif
