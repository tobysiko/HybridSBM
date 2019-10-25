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

#ifndef Cysteine_HH
#define Cysteine_HH
#include "AminoAcid.hh"
#include "../Elements/ATetGroup.hh"
#include "../Elements/DihedralGroup.hh"

namespace prf
{
//! Cystine
    /**
     * Specification of cysteine side chain. SH represented by DihedralGroup.
     * \ingroup aminoacids
     */

    class Cysteine:public AminoAcid
    {

    public:
        Cysteine();
        ~Cysteine();
        void Initialize();

    private:
        ATetGroup g0;
        DihedralGroup g1;
    };
}

#endif
