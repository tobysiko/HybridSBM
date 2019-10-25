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

#include "phenylalanine.hh"

using namespace UnivConstants;

using namespace prf;

Phenylalanine::Phenylalanine() : AminoAcid(F)
{
    node.resize(nnodes=2);
}

Phenylalanine::~Phenylalanine() {}

void Phenylalanine::Initialize()
{
    const double c1=1.52,bXH=1.0;
    double th1,thH;
    th1=113.7*pi/180;
    //double th2=120.0*pi/180;
    thH=108.1*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),atm,mygrp->index(" CD1"));
    g0.SetTheta(thCb,th1,thH,thH);
    g1.SetRootAngle(th1);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetRootLengths(bCaCb,c1);
    node[0]=&g0;
    node[1]=&g1;
    node[0]->AddSubnode(node[1]);
    g0.Initialize();
    g1.Initialize();
    node[0]->SetMobileAtoms(atm[icb+1].UniqueId(),atm[icb+nsatm-1].UniqueId());
    node[1]->SetMobileAtoms(at(" CD1").UniqueId(),atm[icb+nsatm-1].UniqueId());
}
