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

#ifndef BGS_HH
#define BGS_HH
#include "Update.hh"

namespace prf
{
//! Biased Gaussian Steps
    /**
     * Currently this is the most sophisticated update in PROFASI. A set of upto
     * 8 backbone angles are chosen and given a conserted rotation to achieve the
     * following:
     * <ul>
     * <li> A section of the chain somewhere in the middle gets a small deformation
     * </li>
     * <li> The section to one side of the deformed part remains unchanged</li>
     * <li> The section to the other side gets a tiny rigidbody translation and
     * rotation </li>
     * </ul>
     * Therefore, BGS divides the system into three parts, a fixed part, a rigidly
     * moving part and a flexible part, involving atoms around the changed
     * backbone angles.
     *
     * Since this is a semi-local move, the choice of update angles depends on the
     * current configuration of the system. This makes it necessary to modify the
     * Metropolis weight in order to ensure detailed balance. In PROFASI, this is
     * achieved by defining the IntrinsicWeight for all updates. BGS calculates
     * this weight by taking into account the initial and final configurations.
     * Currently all other updates in PROFASI just return 1 for the IntrinsicWeight
     * \ingroup profasi_updates
     */

    class BGS:public prf::Update
    {
    public:
        BGS();
        ~BGS();
        inline void ABGS(double x) { abgs = x; }
        inline void BBGS(double x) { bbgs = x; }
        double intrinsic_weight() const;
        void build_dof_list();
        int perform();
        int revert();
        void print_setup(std::string &st);
    private:
        double abgs, bbgs, bgswt, dph[8];
        int iph[8], nph,dirn;
    };
}

#endif
