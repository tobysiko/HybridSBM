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

#ifndef TrigonalGroup_HH
#define TrigonalGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! A Node with two out-going bonds.
    /**
     * This group is planar, unlike the ATriGroup, which also has two out-going
     * bonds. By planar, we mean, Root, Junction, Branch1 and Branch2 lie in
     * one plane. This is meant to represent groups like CO2, NH2 etc, where
     * the outgoing atoms are identical. The Root-Junction-Branch_i angle is
     * adjustible.
     * \ingroup geometrical_objects
     */

    class TrigonalGroup : public Node
    {
    public:
        TrigonalGroup();
        ~TrigonalGroup();
        void Initialize();
        //! Set the two available theta angles
        inline void SetTheta(double x0, double x1) {
            th0=x0;th1=x1;
        }

        //! Set bond lengths
        inline void SetBondLengths(double x0,double x1,double x2) {
            bprm=x0;bstm=x1;brnc=x2;
        }

        void SetBondLengths(double, double, double, double, double);
        void SetBondAngles(double,double,double,double);
        void SetBranchLength(int i, double x);
        void Create();
    private:
        FindCoord locate_H;
        double bprm,bstm,brnc;
        double th0,th1;
    };
}

#endif
