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

#include "NMethyl.hh"

using namespace UnivConstants;
using std::deque;
using std::pair;

namespace prf
{
    NMethyl::NMethyl() : EndGroup(NME)
    {
        nd = 1;
        na = 0;
    }

    NMethyl::~NMethyl() {}

    void NMethyl::Initialize()
    {
        double thcnh = 121.7 * pi / 180;
        double thccn = 116.6 * pi / 180;
        double bXH=1.0;

        if (bv0 == NULL || bv1 == NULL) {
            cerr <<
            "NMethyl:failed to initialize as the group is not connected\n"
            << "to the C terminal of an Amino Acid\n";
            cerr << "reference 3-vectors are null\n";
            exit(1);
        } else {
            //cout <<"Initialization of NMethyl\n";
            g0.AssignAtoms(rf0, rf1, at(" N  "), at(" H  "),
                           at(" CH3"));
            g1.AssignAtoms(rf1,at(" N  "),at(" CH3"),
                           atm,mygrp->index("1H  "));

            g0.SetBondLengths(AminoAcid::b[2],AminoAcid::b[0],
                              AminoAcid::bNH,AminoAcid::b[1]);
            g0.SetTheta(thccn, thcnh, AminoAcid::theta[0]);
            g0.Initialize();
            g0.LockPhi(0);
            g1.SetBondLengths(AminoAcid::b[0],AminoAcid::b[1],bXH);
            g1.SetRootAngle(AminoAcid::theta[0]);
            g1.Initialize();
            g0.AddSubnode(&g1);
            g0.SetMobileAtoms(at(" H  ").UniqueId(),
                              at("3H  ").UniqueId());
        }
        node.push_back(&g0);
        node.push_back(&g1);
    }

    void NMethyl::Reconstruct()
    {
        atm[0].Pos(rf1.Position() + (*bv1));
        g0.ReCreate();
    }

    void NMethyl::Donor(int i, Dipole & dp)
    {
        dp.SetAtoms(atm[0].UniqueId(), atm[1].UniqueId());
    }

    void NMethyl::Write()
    {
        cout << "NMethyl Group\n";

        for (int i = 0; i < 6; ++i)
            atm[i].Write();
    }

    int NMethyl::rotDof_assign_and_reconstruct(int il, double mgd, int &a0, int &a1)
    {
        g1.AssignPhi(mgd);
        g1.MobileAtoms(a0, a1);
        return 1;
    }

    int NMethyl::rotDof_assign(int il, double mgd)
    {
        g1.RevertPhi(mgd);
        return 1;
    }

    double NMethyl::get_rotDof(int il)
    {
        return g1.Phi();
    }

    int NMethyl::ROTDOF(int il, double mgd, int &a0, int &a1)
    {
        g1.AssignPhi(mgd);
        g1.MobileAtoms(a0, a1);
        return 1;
    }

    int NMethyl::ROTDOFr(int il, double mgd)
    {
        g1.RevertPhi(mgd);
        return 1;
    }

    double NMethyl::ADOF(int il)
    {
        return g1.Phi();
    }

    void NMethyl::ExportConnections(ConnectionsMatrix & aa)
    {
        int egst = atm[0].UniqueId();

        for (int i=0;i<3;++i) {
            aa.set_connection(egst+i,res->Cprime().UniqueId());
            aa.set_connection(egst+i,res->Calpha().UniqueId());
            aa.set_connection(egst+i,res->Cprime().UniqueId()+1);
        }

        g0.ExportConnections(aa);
    }

    void NMethyl::LocPairsatRTdof(int i, deque <pair<int,int> >&lcp)
    {
        switch (i - rtOffset) {

            case 0:
                g1.LocPairs(lcp);
                break;

            default: {
                cout <<
                "asked for local pairs at dof number "
                << (i - rtOffset) << " in N-Methyl\n";
            }
        };
    }

    void NMethyl::BuildConnections()
    {
        g0.set_bases(res->Calpha(),res->Oc());
        g0.BuildConnections();
    }
}
