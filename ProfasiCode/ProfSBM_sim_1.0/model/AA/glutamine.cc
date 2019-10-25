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

#include "glutamine.hh"

using namespace UnivConstants;

using namespace prf;

Glutamine::Glutamine() : AminoAcid(Q)
{
    node.resize(nnodes=4);
}

Glutamine::~Glutamine() {}

void Glutamine::Initialize()
{
    const double c1=1.52,c2=1.52,c3=1.23,c4=1.33,bXH=1.0;
    double th1,th2,th3,th4,thHb,thHg,thHe;
    th1=113.7*pi/180;
    th2=112.6*pi/180;
    th3=121.0*pi/180;
    th4=116.9*pi/180;
    thHb=108.1*pi/180;
    thHg=108.4*pi/180;
    thHe=120.0*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),
                   at(" CD "),at("1HG "),at("2HG "));
    g2.AssignAtoms(Cbeta(),at(" CG "),at(" CD "),
                   at(" OE1"),at(" NE2"));
    g3.AssignAtoms(at(" CG "),at(" CD "),at(" NE2"),
                   at("1HE2"),at("2HE2"));
    g0.SetTheta(thCb,th1,thHb,thHb);
    g1.SetTheta(th1,th2,thHg,thHg);
    g2.SetTheta(th2,th3,th4);
    g3.SetTheta(th4,thHe);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetBondLengths(bCaCb,c1,c2,bXH,bXH);
    g2.SetBondLengths(c1,c2,c3,c4);
    g3.SetBondLengths(c2,c4,bXH);
    g3.LockPhi(0.0);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
    node[3]=&g3;

    for (int j=0;j<3;++j) node[j]->AddSubnode(node[j+1]);

    for (int j=0;j<4;++j) node[j]->Initialize();

    for (int j=0;j<3;++j)
        node[j]->SetMobileAtoms(atm[icb+1+3*j].UniqueId(),
                                at("2HE2").UniqueId());
}


