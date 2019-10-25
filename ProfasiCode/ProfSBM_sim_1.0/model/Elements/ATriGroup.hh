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

#ifndef ATriGroup_HH
#define ATriGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Assymetric Triangular group.
    /**
     * A Node with two outgoing bonds where the incoming and outgoing bonds
     * don't necessarily lie in a plane. Similar to an ATetGroup, with one
     * outgoing bond missing. This is useful while constructing proline. The
     * last CH2 in the side chain, that is connected to the N, is implemented
     * as an ATriGroup
     * \ingroup geometrical_objects
     */

    class ATriGroup : public Node
    {
    public:
        ATriGroup();
        ~ATriGroup();
        void Initialize();
        //! Set independent values for two Branch theta angles
        inline void SetTheta(double x0, double x1,double x2) {
            th0=x0;th1a=x1;th1b=x2;
        }

        //! Set independent values for branch bond lengths
        inline void SetBondLengths(double x0,double x1,double x2,double x3) {
            bprm=x0;bstm=x1;brnca=x2;brncb=x3;
        }

        //! Torsional opening angle between two branches
        inline void SetPhi_ab(double x0) {phba=x0;}

        void Create();
        void SetBondLengths(double, double, double, double, double);
        void SetBranchLength(int i, double x);
        void SetBondAngles(double,double,double,double);
        void SetRelPhi(double,double);
    private:
        FindCoord locate_Ha,locate_Hb;
        double bprm,bstm,brnca,brncb;
        double th0,th1a,th1b,phba;
    };
}

#endif
