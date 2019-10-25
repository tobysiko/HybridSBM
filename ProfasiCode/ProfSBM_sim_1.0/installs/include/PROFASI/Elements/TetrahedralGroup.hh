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

#ifndef TetrahedralGroup_HH
#define TetrahedralGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! The Tetrahedral group
    /**
     * A representation of groups such as the -CH_3. Suppose the group is
     * attached to a certain atom "Root". If Base  is just a hydrogen, all
     * bonds of the carbon atom are identical,  and the angle Root-C-H is
     * arccos(-1/3). But  if Root is a larger atom,  this angle could  be
     * greater. It  will be kept as  an adjustible parameter, "BranchAngle".
     * The three  other atoms attached  to the  carbon are assumed  to be
     * identical.  For the  case where  they are  different, different classes
     * will  be used, like the  ATetGroup class. Another  angle "RootAngle" is
     * required  to  define  a  coordinate   system  with  the  Root  and
     * Axis:  the Base-Root-C angle.
     * \ingroup geometrical_objects
    */

    class TetrahedralGroup : public Node
    {
    public:
        TetrahedralGroup();
        TetrahedralGroup(Atom &,Atom &,Atom &,std::vector<Atom> &,int st);
        ~TetrahedralGroup();
        void Initialize(double,double,double,double,double);
        void Initialize();
        //! Set the Root-Junction-Outgoing angle
        inline void SetBranchAngle(double x) {th12=x;}

        //! Set the Base-Root-Junction angle
        inline void SetRootAngle(double x) {th01=x;}

        //! Set bondlengths Base-Root, Root-Junction, Junction-Branch
        inline void SetBondLengths(double x0, double x1,double x2) {
            bprm=x0;bstm=x1;brnc=x2;
        }

        //! Set bondlength for Junction-Branch, the integer is irrelevant here
        void SetBranchLength(int i, double x);
        void SetBondLengths(double, double, double, double, double);
        void SetBondAngles(double,double,double,double);

        void Create();
    private:
        FindCoord locate_H;
        double bprm,bstm,brnc;
        double th01,th12;
    };
}

#endif
