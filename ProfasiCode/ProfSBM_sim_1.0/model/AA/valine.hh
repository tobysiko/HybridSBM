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

#ifndef Valine_HH
#define Valine_HH
#include "AminoAcid.hh"
#include "../Elements/ATetGroup.hh"
#include "../Elements/TetrahedralGroup.hh"

namespace prf
{
//! Valine
    /**
     * Geometry of valine : an asymmetric tetrahedron (CH), and two tetrahedral
     * nodes with completely symmetric "outgoing" bonds (CH3).
     * \ingroup aminoacids
     */

    class Valine:public AminoAcid
    {

    public:
        Valine();
        ~Valine();
        void Initialize();

    private:
        ATetGroup g0;
        TetrahedralGroup g1, g2;
    };
}

#endif
