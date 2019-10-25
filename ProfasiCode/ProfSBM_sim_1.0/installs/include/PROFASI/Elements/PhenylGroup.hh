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

#ifndef PhenylGroup_HH
#define PhenylGroup_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Represents the phenyl ring with attached hydrogens.
    /**
     * The Phenyl ring with all attached hydrogens is regarded as one rigid
     * block. So, it is easy to calculate positions of all these atoms when
     * the position of two of the carbon atoms in the ring is fixed. It is
     * therefore just one node even if there are many junctions where bonds
     * meet. Root here is the external atom to which the ring is attached.
     * Junction is the Carbon on the ring attached to the Root. There is no
     * branch. The whole ring is regarded as the outgoing part.
     * \ingroup geometrical_objects
     */

    class PhenylGroup : public Node
    {
    public:
        PhenylGroup();
        ~PhenylGroup();
        //! Ininitalize now calculates geometrical constants for phenyl ring
        void Initialize();
        //! Base-Root and Root-Junction lengths
        inline void SetRootLengths(double x0,double x1) {bprm=x0;bstm=x1;}

        //! Base-Root-Junction angle
        inline void SetRootAngle(double x) {th01=x;}

        void Create();
        void ExportConnections(ConnectionsMatrix &aa);
        void LocPairs(std::deque<std::pair<int,int> > & lcp);
        void SetBondLengths(double, double, double, double, double);
        void SetBondAngles(double,double,double,double);

    private:
        FindCoord hexarm;
        double th01;
        double bstm,bprm;
        static double ringbl,bXH;
    };
}

#endif
