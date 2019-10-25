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

#ifndef ATetGroup_HH
#define ATetGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Asymmetric Tetrahedral Group
    /**
     * Mainly represents the CH2 groups in side chains. But it is meant to be
     * any atom connected to 4 others. The three outgoing bonds could have
     * different kinds of atoms, and hence different bondlengths and bond
     * angles.
     * \ingroup geometrical_objects
     */

    class ATetGroup : public Node
    {
    public:
        ATetGroup();
        ATetGroup(Atom &,Atom &,Atom &,std::vector<Atom> &,int st);
        ~ATetGroup();
        //! Assign individually the bond angle for each outgoing bond
        /**
         * The first parameter is the Base-Root-Junction angle, the other 3 are
         * Root-Junction-Branch_i angles
         */
        inline void SetTheta(double t0,double t1,double t2,double t3) {
            th01=t0;th1a=t1;th1b=t2;th1c=t3;
        }

        //! Adjust the locations of the Branches around the incoming bond
        inline void SetPhi(double pab, double pac) {phab=pab;phac=pac;}

        void Initialize();

        void SetBondLengths(double, double, double, double, double);
        //! Set the bondangles base-root-junction and root-junction-branch1,2,3
        void SetBondAngles(double,double,double,double);
        //! Set the phi angle separation between brances 1-2 and 1-3
        void SetRelPhi(double,double);
        //! Set bondlength for branch 0,1 or 2
        void SetBranchLength(int i, double x);
        void Create();
    private:
        FindCoord locate_Ha,locate_Hb,locate_Hc;
        double bprm,bstm,brnca,brncb,brncc;
        double th01;
        double th1a,th1b,th1c;
        double phab,phac;
    };
}

#endif
