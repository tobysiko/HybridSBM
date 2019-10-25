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

#include "isoleucine.hh"

using namespace UnivConstants;

using namespace prf;

Isoleucine::Isoleucine() : AminoAcid(I)
{
    node.resize(nnodes=4);
}

Isoleucine::~Isoleucine() {}

void Isoleucine::Initialize()
{
    const double c1=1.53,c2=1.52,bXH=1.0;
    double th1,th2,thHb,thHg;
    th1=110.4*pi/180;
    th2=113.6*pi/180;
    thHb=109.2*pi/180;
    thHg=108.1*pi/180;
    //double thHd=109.5*pi/180;
    //double thHg2=109.5*pi/180;

    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG1"),at(" HB "),at(" CG2"));
//note reverse order to compensate for the inconsistency between the
//conventional PDB ordering and the chiral ordering of the two gamma
//carbon atms in the side chain.
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG1"),at(" CD1"),at("1HG1"),at("2HG1"));
    g2.AssignAtoms(Calpha(),Cbeta(),at(" CG2"),atm,mygrp->index("1HG2"));
    g3.AssignAtoms(Cbeta(),at(" CG1"),at(" CD1"),atm,mygrp->index("1HD1"));
    g0.AddSubnode(&g1);
    g0.AddSubnode(&g2);
    g1.AddSubnode(&g3);
    g0.SetTheta(thCb,th1,thHb,th1);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,c1);
    g1.SetTheta(th1,th2,thHg,thHg);
    g1.SetBondLengths(bCaCb,c1,c2,bXH,bXH);
    g2.SetRootAngle(th1);
    g2.SetBondLengths(bCaCb,c1,bXH);
    g3.SetRootAngle(th2);
    g3.SetBondLengths(c1,c2,bXH);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
    node[3]=&g3;

    for (int i=0;i<4;++i) node[i]->Initialize();

    g0.SetMobileAtoms(at(" HB ").UniqueId(),at("3HD1").UniqueId());

    g1.SetMobileAtoms(at("1HG1").UniqueId(),at("3HD1").UniqueId());
}
