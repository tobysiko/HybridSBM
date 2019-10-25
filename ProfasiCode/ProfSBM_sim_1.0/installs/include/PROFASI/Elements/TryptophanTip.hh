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

#ifndef TryptophanTip_HH
#define TryptophanTip_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Tip of the side chain of Tryptophan
    /**
     * The two rings, one pentagon the other hexagon of Tryptophan
     */

    class TryptophanTip : public Node
    {
    public:
        TryptophanTip();
        ~TryptophanTip();
        void Initialize();
        void SetBondLengths(double,double,double,double,double);
        void SetBondLengths(double,double,double);
        void Create();
        void ExportConnections(ConnectionsMatrix &aa);
        void LocPairs(std::deque<std::pair<int,int> > & lcp);
    private:
        FindCoord pentarm;
        double th01,th12;
        double bstm,bprm;
        double ringbl,bXH,cu,cv,cw,cw5,cq5,cq6;
    };
}

#endif
