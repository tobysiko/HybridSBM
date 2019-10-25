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

#ifndef Acetyl_HH
#define Acetyl_HH
#include "EndGroup.hh"
#include "../AA/Ligand.hh"
#include "../Elements/ATriGroup.hh"
#include "../Elements/TetrahedralGroup.hh"

namespace prf
{

    class Acetyl:public EndGroup
    {

    public:
        Acetyl();
        ~Acetyl();
        void Initialize();
        void Reconstruct();
        void ExportConnections(ConnectionsMatrix & aa);
        void Acceptor(int i, Dipole & dp);
        void Write();
        void LocPairsatRTdof(int i, std::deque <std::pair<int,int> >&lcp);
        void BuildConnections();
        int rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1);
        int rotDof_assign(int il, double mgd);
        double get_rotDof(int il);
        //OBSOLETE!!
        int ROTDOF(int i, double mgd, int &a0, int &a1);
        int ROTDOFr(int il, double mgd);
        double ADOF(int il);

    private:
        ATriGroup g0;
        TetrahedralGroup g1;
    };
}


#endif
