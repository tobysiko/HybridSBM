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

#ifndef MPhenylGroup_HH
#define MPhenylGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Modified Phenyl Group
    /**
     * Almost identical to the PhenylGroup, but all atoms attached to the
     * ring can now be different. They could be located at different distances.
     * Useful to construct the ring of the Tyrosine. In principle this class
     * could serve to represent PhenylGroup as well. But to keep open this
     * possibility of different atoms is more expensive. So, it should be
     * used only when required.
     * \ingroup geometrical_objects
     */

    class MPhenylGroup : public Node
    {
    public:
        MPhenylGroup();
        ~MPhenylGroup();
        void Initialize();
        inline void SetRootLengths(double x0,double x1) {bprm=x0;bstm=x1;}

        inline void SetRootAngle(double x) {th01=x;}

        void Create();
        void ExportConnections(ConnectionsMatrix &aa);
        void BuildConnections();
        void LocPairs(std::deque<std::pair<int,int> > & lcp);
        void SetBondLengths(double, double, double, double, double);
        void SetBondAngles(double,double,double,double);
        //! Length of the bond connecting the i'th attachment to the ring
        void SetBranchLength(int i, double x);
    private:
        FindCoord hexarm;
        double th01;
        double bstm,bprm;
        double brnc[5];
        static double ringbl,bXH;
    };
}


#endif
