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

#include "cysteine.hh"

using namespace UnivConstants;

using namespace prf;
Cysteine::Cysteine() : AminoAcid(C)
{
    node.resize(nnodes=2);
}

Cysteine::~Cysteine() {}

void Cysteine::Initialize()
{
    const double c1=1.81,bXH=1.0;
    double th,thH,thHg;
    th=113.4*pi/180;
    thH=108.1*pi/180;
    thHg=108.0*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" SG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" SG "),at(" HG "));
    g0.SetTheta(thCb,th,thH,thH);
    g1.SetTheta(th,thHg);
    g0.SetBondLengths(b[1],bCaCb,c1, bXH,bXH);
    g1.SetBondLengths(bCaCb,c1,bXH);
    node[0]=&g0;
    node[1]=&g1;
    node[0]->AddSubnode(node[1]);

    for (int j=0;j<2;++j) node[j]->Initialize();

    node[0]->SetMobileAtoms(at("1HB ").UniqueId(),at(" HG ").UniqueId());
}

