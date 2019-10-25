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

#include "histidine.hh"

using namespace UnivConstants;

using namespace prf;

Histidine::Histidine() : AminoAcid(H)
{
    node.resize(nnodes=2);
}

Histidine::~Histidine() {}

void Histidine::Initialize()
{
    const double c1=1.52,c2=1.35,bXH=1.0;
    double th1,thH;
    th1=113.2*pi/180;
    //double th2=126.0*pi/180; //now taken care of in HistidineTip class
    thH=108.2*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),atm,mygrp->index(" CD2"));
    g0.SetTheta(thCb,th1,thH,thH);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetBondLengths(bCaCb,c1,c2);
    node[0]=&g0;
    node[1]=&g1;
    node[0]->AddSubnode(node[1]);

    for (int j=0;j<2;++j) node[j]->Initialize();

    node[0]->SetMobileAtoms(atm[icb+1].UniqueId(),
                            atm[icb+nsatm-1].UniqueId());
}
