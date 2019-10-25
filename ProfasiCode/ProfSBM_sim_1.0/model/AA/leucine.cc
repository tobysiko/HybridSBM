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

#include "leucine.hh"

using namespace UnivConstants;

using namespace prf;

Leucine::Leucine() : AminoAcid(L) {node.resize(nnodes=4);}

Leucine::~Leucine() {}

void Leucine::Initialize()
{
    const double c1=1.53,c2=1.52,bXH=1.0;
    double th1,th2,thHb,thHg;
    th1=117.1*pi/180;
    th2=110.1*pi/180;
    thHb=107.0*pi/180;
    thHg=109.3*pi/180;
    //double thHd=109.5*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),
                   at(" CD1"),at(" CD2"),at(" HG "));
    g2.AssignAtoms(Cbeta(),at(" CG "),at(" CD1"),atm,mygrp->index("1HD1"));
    g3.AssignAtoms(Cbeta(),at(" CG "),at(" CD2"),atm,mygrp->index("1HD2"));
    g0.AddSubnode(&g1);
    g1.AddSubnode(&g2);
    g1.AddSubnode(&g3);
    g0.SetTheta(thCb,th1,thHb,thHb);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetTheta(th1,th2,th2,thHg);
    g1.SetBondLengths(bCaCb,c1,c2,c2,bXH);
    g2.SetRootAngle(th2);
    g2.SetBondLengths(c1,c2,bXH);
    g3.SetRootAngle(th2);
    g3.SetBondLengths(c1,c2,bXH);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
    node[3]=&g3;

    for (int i=0;i<4;++i) node[i]->Initialize();

    g1.SetMobileAtoms(at(" HG ").UniqueId(),at("3HD2").UniqueId());

    g0.SetMobileAtoms(at("1HB ").UniqueId(),at("3HD2").UniqueId());
}
