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

#include "Acetyl.hh"

using namespace UnivConstants;
using std::deque;
using std::pair;

namespace prf
{
    Acetyl::Acetyl() : EndGroup(ACE)
    {
        nd = 0;
        na = 1;
    }

    Acetyl::~Acetyl() {}

    void Acetyl::Initialize()
    {
        double thnco = 121.7 * pi / 180;
        double thncc = 116.6 * pi / 180;
        double bXH=1.0;

        if (bv0 == NULL || bv1 == NULL) {
            cerr <<
            "Acetyl:failed to initialize as the group is not connected\n"
            << "to the N terminal of an Amino Acid\n";
            cerr << "reference 3-vectors are null\n";
            exit(1);
        } else {
            //cout <<"Initialization of Acetyl\n";
            g0.AssignAtoms(rf0, rf1, at(" C  "), at(" O  "),
                           at(" CH3"));
            g1.AssignAtoms(rf1,at(" C  "),at(" CH3"),
                           atm,mygrp->index("1H  "));

            g0.SetBondLengths(AminoAcid::b[1],AminoAcid::b[0],
                              AminoAcid::bCO,AminoAcid::b[2]);
            g0.SetTheta(AminoAcid::theta[0], thnco, thncc);
            g0.Initialize();
            g0.LockPhi(0);
            g1.SetBondLengths(AminoAcid::b[0],AminoAcid::b[2],bXH);
            g1.SetRootAngle(thncc);
            g1.Initialize();
            g0.AddSubnode(&g1);
            g0.SetMobileAtoms(at(" O  ").UniqueId(),
                              at("3H  ").UniqueId());
        }
        node.push_back(&g0);
        node.push_back(&g1);
    }

    void Acetyl::Reconstruct()
    {
        atm[0].Pos(rf1.Position() - (*bv1));
        g0.ReCreate();
    }

    void Acetyl::Acceptor(int i, Dipole & dp)
    {
        dp.SetAtoms(atm[1].UniqueId(), atm[0].UniqueId());
    }

    void Acetyl::Write()
    {
        cout << "Acetyl Group\n";

        for (int i = 0; i < 6; ++i)
            atm[i].Write();
    }

    int Acetyl::rotDof_assign_and_reconstruct(int il, double mgd, int &a0, int &a1)
    {
        g1.AssignPhi(mgd);
        g1.MobileAtoms(a0, a1);
        return 1;
    }

    int Acetyl::rotDof_assign(int il, double mgd)
    {
        g1.RevertPhi(mgd);
        return 1;
    }

    double Acetyl::get_rotDof(int il)
    {
        return g1.Phi();
    }

    int Acetyl::ROTDOF(int il, double mgd, int &a0, int &a1)
    {
        g1.AssignPhi(mgd);
        g1.MobileAtoms(a0, a1);
        return 1;
    }

    int Acetyl::ROTDOFr(int il, double mgd)
    {
        g1.RevertPhi(mgd);
        return 1;
    }

    double Acetyl::ADOF(int il)
    {
        return g1.Phi();
    }

    void Acetyl::ExportConnections(ConnectionsMatrix & aa)
    {
        int aast = res->Nitrogen().UniqueId();
        int egst = atm[0].UniqueId();

        for (int i = 0; i < 3; ++i) {
            if (res->OLC() != P&& res->OLC() !=DPR) {
                for (int j = 0; j < 3; ++j)
                    aa.set_connection(egst + i,aast + j);
            } else {
                for (int j = 0; j < 13; ++j)
                    aa.set_connection(egst + i,aast + j);
            }
        }

        if (res->OLC() != P && res->OLC() !=DPR) {
            aa.set_connection(egst, res->Hca().UniqueId());
            aa.set_connection(egst, res->Cbeta().UniqueId());
            aa.set_connection(egst,res->Cprime().UniqueId());
        }

        g0.ExportConnections(aa);
    }

    void Acetyl::LocPairsatRTdof(int i, deque <pair<int,int> >&lcp)
    {
        switch (i - rtOffset) {

            case 0:
                g1.LocPairs(lcp);
                break;

            default: {
                cout <<
                "asked for local pairs at dof number "
                << (i - rtOffset) << " in acetyl\n";
            }
        };
    }

    void Acetyl::BuildConnections()
    {
        if (res->OLC() != P && res->OLC() !=DPR) {
            g0.set_bases(res->Calpha(),res->at(" H  "));
        } else {
            g0.set_bases(res->Calpha(),res->at(" CD "));
        }

        g0.BuildConnections();
    }
}
