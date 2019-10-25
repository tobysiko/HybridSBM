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

#ifndef Succinyl_HH
#define Succinyl_HH
#include "EndGroup.hh"
#include "../AA/Ligand.hh"
#include "../Elements/ATriGroup.hh"
#include "../Elements/ATetGroup.hh"
#include "../Elements/TrigonalGroup.hh"

namespace prf
{

    class Succinyl:public EndGroup
    {

    public:
        Succinyl();
        ~Succinyl();
        void Allocate();
        void Initialize();
        void Reconstruct();
        void ExportConnections(ConnectionsMatrix & aa);
        int ROTDOF(int i, double mgd, int &a0, int &a1);
        int ROTDOFr(int il, double mgd);
        double ADOF(int il);
        void Acceptor(int i, Dipole & dp);
        void Write();
        void LocPairsatRTdof(int i, std::deque<std::pair<int,int> >&lcp);
        void BuildConnections();
        int rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1);
        int rotDof_assign(int il, double mgd);
        double get_rotDof(int il);

    private:
        ATriGroup g0;
        ATetGroup g1, g2;
        TrigonalGroup g3;
    };
}


#endif
