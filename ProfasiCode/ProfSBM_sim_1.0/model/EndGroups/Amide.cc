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

#include "Amide.hh"

namespace prf
{
    Amide::Amide() : EndGroup(NH2)
    {
        nd = 2;
        na = 0;
    }

    Amide::~Amide() {}

    void Amide::Initialize()
    {
        if (bv0 == NULL || bv1 == NULL) {
            prf::cerr <<"Amide: failed to initialize as the "
            <<"Amide group is not connected\n"
            << "to the C terminal of an Amino Acid\n"
            <<"reference 3-vectors are null\n";
            exit(1);
        } else {
            //cout <<"Initialization of Amide\n";
            grp.AssignAtoms(rf0, rf1, atm[0], atm[1], atm[2]);
            grp.SetBondLengths(AminoAcid::b[2],
                               AminoAcid::b[0],
                               AminoAcid::bNH);
            grp.SetTheta(AminoAcid::theta[2], UnivConstants::twoPid3);
            grp.Initialize();
            grp.LockPhi(0);
        }
        node.push_back(&grp);
    }

    void Amide::Reconstruct()
    {
        atm[0].Pos(rf1.Position() + (*bv1));
        grp.Create();
    }

    void Amide::Donor(int i, Dipole & dp)
    {
        i = i % nd;

        if (i == 0)
            dp.SetAtoms(atm[0].UniqueId(), atm[1].UniqueId());
        else
            dp.SetAtoms(atm[0].UniqueId(), atm[2].UniqueId());
    }

    void Amide::Write()
    {
        prf::cout << "Amide Group\n";

        for (int i=0;i<3;++i) atm[i].Write();
    }

    int Amide::ROTDOF(int il, double mgd, int &a0, int &a1)
    {
        if (il == rtOffset)
            grp.AssignPhi(mgd);

        grp.MobileAtoms(a0, a1);

        return 1;
    }

    int Amide::ROTDOFr(int il, double mgd)
    {
        if (il == rtOffset)
            grp.RevertPhi(mgd);

        return 1;
    }

    double Amide::ADOF(int il)
    {
        return grp.Phi();
    }

    void Amide::BuildConnections()
    {
        grp.set_bases(res->Calpha(),res->Oc());
    }

    void Amide::ExportConnections(ConnectionsMatrix & aa)
    {
        grp.ExportConnections(aa);
        aa.set_connection(res->Nitrogen().UniqueId(),atm[0].UniqueId());
        aa.set_connection(res->Hca().UniqueId(),atm[0].UniqueId());
        aa.set_connection(res->Cbeta().UniqueId(),atm[0].UniqueId());
    }
}
