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

#include "valine.hh"

using namespace UnivConstants;

using namespace prf;

Valine::Valine() : AminoAcid(V) {node.resize(3);nnodes=3;}

Valine::~Valine() {}

void Valine::Initialize()
{
    const double c1=1.52,bXH=1.0;
    double th,thH;
    th=110.7*pi/180;
    thH=109.1*pi/180;
    //double thH2=109.5*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),at(" CG1"),at(" CG2"),at(" HB "));
    g0.SetTheta(thCb,th,th,thH);
    g0.SetBondLengths(b[1],bCaCb,c1,c1,bXH);
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG1"),atm,mygrp->index("1HG1"));
    g1.SetRootAngle(th);
    g1.SetBondLengths(bCaCb,c1,bXH);
    g2.AssignAtoms(Calpha(),Cbeta(),at(" CG2"),atm,mygrp->index("1HG2"));
    g2.SetRootAngle(th);
    g2.SetBondLengths(bCaCb,c1,bXH);
    g0.Initialize();
    g1.Initialize();
    g2.Initialize();
    g0.AddSubnode(&g1);
    g0.AddSubnode(&g2);
    g0.SetMobileAtoms(at(" HB ").UniqueId(),at("3HG2").UniqueId());
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
}
