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

#ifndef Proline_HH
#define Proline_HH
#include "AminoAcid.hh"
#include "../Elements/ATetGroup.hh"
#include "../Elements/ATriGroup.hh"

namespace prf
{
//! Proline
    /**
     * Specification of proline side chain. ATriGroup is a trigonal node with two
     * out-going bonds, but the incomming and out-going bonds neeed not lie in a
     * plane. The last CH2 of the proline side-chain
     * which is covalently bonded to the N is represented by this group. The
     * model ignores the puckered nature of the proline side chain, so that all
     * three carbon atoms are assumed to lie in a plane. The side chain is
     * constructed starting from the Cbeta as usual. Since we assume all proline
     * side chain chi angles as well as the proline backbone phi angle to be fixed
     * the Cdelta always appears at a fixed distance from the N.
     * \ingroup aminoacids
     */

    class Proline:public AminoAcid
    {

    public:
        Proline();
        Proline(ChiralityType ch);
        ~Proline();
        void Initialize();
        void Reconstruct();

    private:
        ATetGroup g0, g1;
        ATriGroup g2;
    };
}

#endif
