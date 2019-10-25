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

#ifndef Pivot_HH
#define Pivot_HH
#include "Update.hh"

namespace prf
{
//! Pivot move twists the protein about one backbone degree of freedom
    /**
     * Only one Ramachandran phi or psi angle is chosen and changed. This is
     * equivalent to a rigid body rotation to a part of the system. Implementation
     * in this way is much faster than reconstructing the whole chain from the
     * internal coordinates after changing one.
     *
     * Pivot divides the system into two parts that rotate rigidly relative to
     * each other.
     *
     * When a backbone angle is changed, either the part of the molecule from the
     * that point towards the C-terminus or from that point to the N-terminus
     * could be rotated with respect to the rest. This defines a "direction" for the
     * update. The first of the above mentioned direction is defined to have
     * value 0 the second 1. The direction is chosen so that the shorter chain
     * segment is rotated with greater probability.
     * \ingroup profasi_updates
     *
     */

    class Pivot:public Update
    {
    public:
        Pivot();
        ~Pivot();
        int perform();
        void build_dof_list();
        int revert();
    private:
        int dirn;
        double scl;
    };
}

#endif
