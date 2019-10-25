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

#include "asparagine.hh"

using namespace UnivConstants;

using namespace prf;

Asparagine::Asparagine() : AminoAcid(N)
{
    node.resize(nnodes=3);
}

Asparagine::~Asparagine() {}

void Asparagine::Initialize()
{
    const double c1=1.52,c2=1.23,c3=1.33,bXH=1.0;
    double th1,th2,th3,thH,thHd;
    th1=112.6*pi/180;
    th2=120.9*pi/180;
    th3=117.0*pi/180;
    thH=108.4*pi/180;
    thHd=120.0*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),
                   at(" OD1"),at(" ND2"));
    g2.AssignAtoms(Cbeta(),at(" CG "),at(" ND2"),
                   at("1HD2"),at("2HD2"));
    g0.SetTheta(thCb,th1,thH,thH);
    g1.SetTheta(th1,th2,th3);
    g2.SetTheta(th3,thHd);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetBondLengths(bCaCb,c1,c2,c3);
    g2.SetBondLengths(c1,c3,bXH);
    g2.LockPhi(0.0);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
    node[0]->AddSubnode(node[1]);
    node[1]->AddSubnode(node[2]);

    for (int j=0;j<3;++j) node[j]->Initialize();

    g0.SetMobileAtoms(atm[icb+1].UniqueId(),at("2HD2").UniqueId());

    g1.SetMobileAtoms(at(" OD1").UniqueId(),at("2HD2").UniqueId());
}


