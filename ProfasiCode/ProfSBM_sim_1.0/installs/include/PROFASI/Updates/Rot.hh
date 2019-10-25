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

#ifndef Rot_HH
#define Rot_HH
#include "Update.hh"

namespace prf
{
//! Rotation of a single side chain torsional degree of freedom
    /**
     * The simplest update in PROFASI. Side-chains of all kinds of residues can
     * be regarded as short compared to typical backbone lengths. So, after
     * changing the side chain angle, the side chain is reconstructed down-stream
     * from the values of the chi angles.
     * \ingroup profasi_updates
     */

    class Rot:public Update
    {
    public:
        Rot();
        ~Rot();
        int perform();
        void build_dof_list();
    };
}

#endif
