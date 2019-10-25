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

#ifndef Tyrosine_HH
#define Tyrosine_HH
#include "AminoAcid.hh"
#include "../Elements/ATetGroup.hh"
#include "../Elements/MPhenylGroup.hh"
#include "../Elements/DihedralGroup.hh"

namespace prf
{
//! Tyrosine
    /**
     * Specification of tyrosine side chain. Here we use an MPhenylGroup which
     * creates a phenyl ring, but the atoms attached to the ring carbon atoms
     * could be connected with bonds of different length. This feature is used
     * for the OH group attached to the ring.
     * \ingroup aminoacids
     */

    class Tyrosine:public AminoAcid
    {

    public:
        Tyrosine();
        ~Tyrosine();
        void Initialize();

    private:
        ATetGroup g0;
        MPhenylGroup g1;
        DihedralGroup g2;
    };
}

#endif
