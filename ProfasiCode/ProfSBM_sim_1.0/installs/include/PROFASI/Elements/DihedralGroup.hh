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

#ifndef DihedralGroup_HH
#define DihedralGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! The dihedral group or boomerang
    /**
     * This group with only one outgoing bond only occurs in Methionine, Serine
     * Threonine and Tyrosine
     * \ingroup geometrical_objects
     */

    class DihedralGroup : public Node
    {
    public:
        DihedralGroup();
        ~DihedralGroup();
        void Initialize();
        //! Base-Root-Junction, and the only other theta angle in this case
        inline void SetTheta(double x0, double x1) {
            th0=x0;th1=x1;
        }

        //! All bond lengths
        inline void SetBondLengths(double x0,double x1,double x2) {
            bprm=x0;bstm=x1;brnc=x2;
        }

        void SetBondLengths(double, double, double, double, double);
        void SetBranchLength(int i, double x);
        void SetBondAngles(double,double,double,double);
        void Create();
    private:
        FindCoord locate_H;
        double bprm,bstm,brnc;
        double th0,th1;
    };
}

#endif
