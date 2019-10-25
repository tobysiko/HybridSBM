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

#ifndef Translation_HH
#define Translation_HH
#include "Update.hh"

namespace prf
{
//! Rigid body translation of a whole chain
    /**
     * Simple block translation of a whole chain by a random vector. The size of
     * the translation vector is kept small compared to the box size.
     * \ingroup profasi_updates
     */

    class Translation:public Update
    {
    public:
        Translation();
        ~Translation();
        int perform();
        void build_dof_list();
        void setScale(double gscl);
    private:
        double scl;
    };
}

#endif
