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

#ifndef Amide_HH
#define Amide_HH
#include "EndGroup.hh"
#include "../AA/Ligand.hh"
#include "../Elements/TrigonalGroup.hh"

namespace prf
{

    class Amide:public EndGroup
    {

    public:
        Amide();
        ~Amide();
        void Initialize();
        void Reconstruct();
        void ExportConnections(ConnectionsMatrix & aa);
        int ROTDOF(int i, double mgd, int &a0, int &a1);
        int ROTDOFr(int il, double mgd);
        double ADOF(int il);
        void Donor(int i, Dipole & dp);
        void Write();
        void BuildConnections();

    private:
        TrigonalGroup grp;
    };
}


#endif
