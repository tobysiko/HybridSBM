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

#include "Succinyl.hh"

using namespace UnivConstants;
using std::deque;
using std::pair;

namespace prf
{
    Succinyl::Succinyl() : EndGroup(SUC)
    {
        nd = 0;
        na = 1;
        //nd and na refer to the "backbone" dipoles. The two CO
        //groups at the end of Succinyl are treated as the end of
        //the side chain of glutamic acid.
    }

    void Succinyl::Allocate()
    {
        EndGroup::Allocate();

        for (int i = 0; i < natms; ++i)
            atm[i].UniqueId(natms - 1 - i);
    }

    Succinyl::~Succinyl() {}

    void Succinyl::Initialize()
    {
        double thnco = 121.7 * pi / 180, thHg =
                           108.0 * pi / 180, thHd = 108.2 * pi / 180;
        double thncc =
            116.6 * pi / 180, th2(113.2 * pi / 180),
            th3(118.5 * pi / 180);
        double b1(1.52), b2(1.25), bXH(1.0);

        if (bv0 == NULL || bv1 == NULL) {
            prf::cerr <<
            "Succinyl:failed to initialize as the group is not connected\n"
            << "to the N terminal of an Amino Acid\n";
            prf::cerr << "reference 3-vectors are null\n";
            exit(1);
        } else {
            //prf::cout <<"Initialization of Succinyl\n";
            g0.AssignAtoms(rf0, rf1, at(" CO "),
                           at(" O  "),at(" C3 "));
            g1.AssignAtoms(rf1,at(" CO "),at(" C3 "),
                           atm,mygrp->index("1H3 "));
            g2.AssignAtoms(at(" CO "),at(" C3 "),at(" C2 "),
                           atm, mygrp->index("1H2 "));
            g3.AssignAtoms(at(" C3 "),at(" C2 "),at(" CO2"),
                           at(" O11"),at(" O12"));
            g0.SetBondLengths(AminoAcid::b[1],AminoAcid::b[0],
                              AminoAcid::bCO, b1);
            g1.SetBondLengths(AminoAcid::b[0], b1, bXH, bXH,b1);
            g2.SetBondLengths(b1, b1, bXH, bXH, b1);
            g3.SetBondLengths(b1, b1, b2);
            g0.SetTheta(AminoAcid::theta[0], thnco, thncc);
            g1.SetTheta(thncc, thHg, thHg, th2);
            g2.SetTheta(th2, thHd, thHd, th2);
            g3.SetTheta(th2, th3);
            g0.Initialize();
            g0.LockPhi(0);
            g1.Initialize();
            g2.Initialize();
            g3.Initialize();
            g0.AddSubnode(&g1);
            g1.AddSubnode(&g2);
            g2.AddSubnode(&g3);
            g0.SetMobileAtoms(at(" O  ").UniqueId(),
                              at(" O12").UniqueId());
            g1.SetMobileAtoms(at("1H3 ").UniqueId(),
                              at(" O12").UniqueId());
            g2.SetMobileAtoms(at("1H2 ").UniqueId(),
                              at(" O12").UniqueId());
        }
        node.push_back(&g0);
        node.push_back(&g1);
        node.push_back(&g2);
        node.push_back(&g3);
    }

    void Succinyl::Reconstruct()
    {
        atm[0].Pos(rf1.Position() - (*bv1));
        g0.ReCreate();
    }

    void Succinyl::Acceptor(int i, Dipole & dp)
    {
        dp.SetAtoms(at(" O  ").UniqueId(),at(" CO ").UniqueId());
    }

    void Succinyl::Write()
    {
        prf::cout << "Succinyl Group\n";

        for (int i = 0; i < natms; ++i)
            atm[i].Write();
    }

    int Succinyl::rotDof_assign_and_reconstruct(int il, double mgd, int &a0, int &a1)
    {
        switch (il % 3) {

            case 0:
                g1.AssignPhi(mgd);
                g1.MobileAtoms(a0, a1);
                break;

            case 1:
                g2.AssignPhi(mgd);
                g2.MobileAtoms(a0, a1);
                break;

            case 2:

            default:
                g3.AssignPhi(mgd);
                g3.MobileAtoms(a0, a1);
                break;
        };

        return 1;
    }

    int Succinyl::rotDof_assign(int il, double mgd)
    {
        switch (il % 3) {

            case 0:
                g1.RevertPhi(mgd);
                return 1;

            case 1:
                g2.RevertPhi(mgd);
                return 1;

            case 2:

            default:
                g3.RevertPhi(mgd);
                return 1;
        };
    }

    double Succinyl::get_rotDof(int il)
    {
        switch (il % 3) {

            case 0:
                return g1.Phi();

            case 1:
                return g2.Phi();

            case 2:

            default:
                return g3.Phi();
        };
    }

    int Succinyl::ROTDOF(int il, double mgd, int &a0, int &a1)
    {
        return rotDof_assign_and_reconstruct(il,mgd,a0,a1);
    }

    int Succinyl::ROTDOFr(int il, double mgd)
    {
        return rotDof_assign(il,mgd);
    }

    double Succinyl::ADOF(int i1)
    {
        return get_rotDof(i1);
    }

    void Succinyl::ExportConnections(ConnectionsMatrix & aa)
    {
        int aast = res->Nitrogen().UniqueId();
        int egst = atm[0].UniqueId();

        for (int i = 0; i < 3; ++i) {
            if (res->OLC() != P && res->OLC() !=DPR) {
                for (int j = 0; j < 3; ++j)
                    aa.set_connection(egst + i,aast + j);
            } else {
                for (int j = 0; j < 13; ++j)
                    aa.set_connection(egst + i,aast + j);
            }
        }

        if (res->OLC() != P&& res->OLC() !=DPR) {
            aa.set_connection(egst, res->Hca().UniqueId());
            aa.set_connection(egst, res->Cbeta().UniqueId());
            aa.set_connection(egst,res->Cprime().UniqueId());
        }

        g0.ExportConnections(aa);
    }

    void Succinyl::LocPairsatRTdof(int i, deque<pair<int,int> >&lcp)
    {
        switch (i - rtOffset) {

            case 0:
                g1.LocPairs(lcp);
                break;

            case 1:
                g2.LocPairs(lcp);
                break;

            case 2:
                g3.LocPairs(lcp);
                break;

            default: {
                prf::cout <<
                "asked for local pairs at dof number "
                << (i - rtOffset) << " in succinyl\n";
            }
        };
    }

    void Succinyl::BuildConnections()
    {
        if (res->OLC() != P&& res->OLC() !=DPR) {
            g0.set_bases(res->Calpha(),res->at(" H  "));
        } else {
            g0.set_bases(res->Calpha(),res->at(" CD "));
        }

        g0.BuildConnections();
    }
}
