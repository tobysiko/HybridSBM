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

#include "proline.hh"

using namespace UnivConstants;

using namespace prf;

Proline::Proline() : AminoAcid(P)
{
    node.resize(nnodes=3);
    Ca_Chirality(LEV);
}

Proline::Proline(ChiralityType ch) : AminoAcid(ch==LEV?P:DPR)
{
    node.resize(nnodes=3);
    Ca_Chirality(ch);
}

Proline::~Proline() {}

void Proline::Initialize()
{
    const double c1=1.51,c2=1.51,bXH=1.0;
    double th1,th2,thHb,thHg,thHd,chi0,chi1,chi2;
    th1=103.3*pi/180;
    th2=110.8*pi/180;
    thHb=111.6*pi/180;
    thHg=109.0*pi/180;
    thHd=110.7*pi/180;
    chi0=2*pi/3;
    chi1=2*pi/3;
    chi2=2*pi/3;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at("1HB "),at("2HB "),at(" CG "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),
                   at("1HG "),at("2HG "),at(" CD "));
    g2.AssignAtoms(Cbeta(),at(" CG "),at(" CD "),
                   at("1HD "),at("2HD "));
    g0.SetTheta(thCb,thHb,thHb,th1);
    g1.SetTheta(th1,thHg,thHg,th2);
    g2.SetTheta(th2,thHd,thHd);
    g0.LockPhi(chi0);
    g1.LockPhi(chi1);
    g2.LockPhi(chi2);
    g2.SetPhi_ab(120.0*pi/180);
    g0.SetBondLengths(b[1],bCaCb,bXH,bXH,c1);
    g1.SetBondLengths(bCaCb,c1,bXH,bXH,c2);
    g2.SetBondLengths(c1,c2,bXH,bXH);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;

    for (int j=0;j<2;++j) node[j]->AddSubnode(node[j+1]);

    for (int j=0;j<3;++j) node[j]->Initialize();

    for (int j=0;j<3;++j) node[j]->SetMobileAtoms(atm[ibbca+3+3*j].UniqueId(),
                atm[ibbca+10].UniqueId());
}

void Proline::Reconstruct()
{
    if (hasNTerminal()) {
        atm[natms-2].Pos(Nitrogen().Pos() +
            locate_H2n(theBB->bond(BBloc), theBB->bond(BBloc+1)));
        atm[natms-1].Pos(Nitrogen().Pos() +
            locate_H3n(theBB->bond(BBloc), theBB->bond(BBloc+1)));
    }

    if (chirality==LEV) {   /* Normal L Proline */
        atm[ibbca+1].Pos(Calpha().Pos() +
            locate_Ha(theBB->bond(BBloc+2), theBB->bond(BBloc+1)));
        atm[ibbca+2].Pos(Calpha().Pos() +
            locate_Cb(theBB->bond(BBloc+2), theBB->bond(BBloc+1)));
    } else { /* D Proline */
        atm[ibbca+1].Pos(Calpha().Pos() +
            locate_dHa(theBB->bond(BBloc+2), theBB->bond(BBloc+1)));
        atm[ibbca+2].Pos(Calpha().Pos() +
            locate_dCb(theBB->bond(BBloc+2), theBB->bond(BBloc+1)));
    }

    atm[ibbc+1].Pos(Cprime().Pos() +
            locate_Oc(theBB->bond(BBloc+2), theBB->bond(BBloc+3)));

    if (hasCTerminal()) {
        atm[ibbc+2].Pos(Cprime().Pos() +
            locate_O2c(theBB->bond(BBloc+2), theBB->bond(BBloc+3)));
    }

    g0.ReCreate();
}

