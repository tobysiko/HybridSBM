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

#ifndef HistidineTip_HH
#define HistidineTip_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Pentagon at the tip of the side chain of Histidine
    /**
     * Exclusively meant to be used in Histidine. Places the Pentagon in space.
     * \ingroup geometrical_objects
     */

    class HistidineTip : public Node
    {
    public:
        HistidineTip();
        ~HistidineTip();
        void Initialize();
        void SetBondLengths(double, double, double,double,double);
        void SetBondLengths(double, double, double);
        void Create();
        void ExportConnections(ConnectionsMatrix &aa);
        void LocPairs(std::deque<std::pair<int,int> > & lcp);
    private:
        FindCoord pentarm;
        double th01,th12;
        double bstm,bprm;
        double ringbl,bXH,cu,cv,cw,cq;
    };
}


#endif
