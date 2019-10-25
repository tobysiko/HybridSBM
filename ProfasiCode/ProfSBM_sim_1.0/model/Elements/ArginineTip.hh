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

#ifndef ArginineTip_HH
#define ArginineTip_HH
#include "Node.hh"
#include "../Aux/FindCoord.hh"
#include "../Aux/Constants.hh"

namespace prf
{
    //! Representation for the rigid tip of the Arginine side chain
    /**
     * Tip of the side chain of Arginine, exclusively for Arginine, so,
     *  the angles and bond lengths are not adjustible.
     * \ingroup geometrical_objects
     */

    class ArginineTip : public Node
    {
    public:
        ArginineTip();
        ArginineTip(Atom &,Atom &,Atom &,std::vector<Atom> &,int st);
        ~ArginineTip();

        void Initialize();
        void Create();
        void LocPairs(std::deque<std::pair<int,int> > & lcp);
        void ExportConnections(ConnectionsMatrix &aa);

    private:
        FindCoord locate_C,locate_N;
        double bprm,bstm,brncnc,brnccn,brncnh;
        double th01;
        double th12,th23,th34;
        double nhdcn,nhdnc;
    };
}

#endif
