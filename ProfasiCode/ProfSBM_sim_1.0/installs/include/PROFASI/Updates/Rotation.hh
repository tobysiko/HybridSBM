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

#ifndef Rotation_HH
#define Rotation_HH
#include "Update.hh"

namespace prf
{
//! Rigid body rotation of a whole chain
    /**
     * Rotates the whole chain by a random angle, about its center of mass.
     * Normally this would not change the energy in a single chain system, unless
     * there is some external field. We have used it mostly for multiple-chain
     * systems.
     * \ingroup profasi_updates
     *
     */

    class Rotation:public Update
    {
    public:
        Rotation();
        ~Rotation();
        int perform();
        void build_dof_list();
        void setScale(double gscl);
    private:
        double scl;
    };
}

#endif
