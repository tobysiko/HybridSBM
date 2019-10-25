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

#include "methionine.hh"

using namespace UnivConstants;

using namespace prf;

Methionine::Methionine() : AminoAcid(M)
{
    node.resize(nnodes=4);
}

Methionine::~Methionine() {}

void Methionine::Initialize()
{
    const double c1=1.52,c2=1.81,c3=1.79,bXH=1.0;
    double th1,th2,th3,thHb,thHg;
    th1=113.5*pi/180;
    th2=111.9*pi/180;
    th3=100.5*pi/180;
    thHb=108.1*pi/180;
    thHg=108.7*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),
                   at(" SD "),at("1HG "),at("2HG "));
    g2.AssignAtoms(Cbeta(),at(" CG "),at(" SD "),at(" CE "));
    g3.AssignAtoms(at(" CG "),at(" SD "),at(" CE "),atm,mygrp->index("1HE "));

    g0.SetTheta(thCb,th1,thHb,thHb);
    g1.SetTheta(th1,th2,thHg,thHg);
    g2.SetTheta(th2,th3);
    g3.SetRootAngle(th3);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetBondLengths(bCaCb,c1,c2,bXH,bXH);
    g2.SetBondLengths(c1,c2,c3);
    g3.SetBondLengths(c2,c3,bXH);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
    node[3]=&g3;

    for (int j=0;j<3;++j) node[j]->AddSubnode(node[j+1]);

    for (int j=0;j<4;++j) node[j]->Initialize();

    for (int j=0;j<2;++j) node[j]->SetMobileAtoms(atm[ibbca+3+3*j].UniqueId(),atm[ibbca+12].UniqueId());

    node[2]->SetMobileAtoms(atm[ibbca+9].UniqueId(),atm[ibbca+12].UniqueId());
}

