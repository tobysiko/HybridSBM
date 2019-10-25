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

#include "arginine.hh"

using namespace UnivConstants;

using namespace prf;

Arginine::Arginine() : AminoAcid(R) {node.resize(nnodes=4);}

Arginine::~Arginine() {}

void Arginine::Initialize()
{
    const double c1=1.52,c2=1.52,c3=1.46,bXH=1.0;
    double th1,th2,th3,thHb,thHg,thHd;
    th1=113.7*pi/180;
    th2=th3=111.5*pi/180;
    thHb=108.1*pi/180;
    thHg=thHd=108.8*pi/180;
    g0.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),
                   at(" CG "),at("1HB "),at("2HB "));
    g1.AssignAtoms(Calpha(),Cbeta(),at(" CG "),
                   at(" CD "),at("1HG "),at("2HG "));
    g2.AssignAtoms(Cbeta(),at(" CG "),at(" CD "),
                   at(" NE "),at("1HD "),at("2HD "));
    g3.AssignAtoms(at(" CG "),at(" CD "),at(" NE "),
                   atm,mygrp->index(" HE "));
    g0.SetTheta(thCb,th1,thHb,thHb);
    g1.SetTheta(th1,th2,thHg,thHg);
    g2.SetTheta(th2,th3,thHd,thHd);
    g0.SetBondLengths(b[1],bCaCb,c1,bXH,bXH);
    g1.SetBondLengths(bCaCb,c1,c2 ,bXH,bXH);
    g2.SetBondLengths(c1,c2,c3,bXH,bXH);
    node[0]=&g0;
    node[1]=&g1;
    node[2]=&g2;
    node[3]=&g3;

    for (int j=0;j<3;++j) node[j]->AddSubnode(node[j+1]);

    for (int j=0;j<4;++j) node[j]->Initialize();

    for (int j=0;j<3;++j)
        node[j]->SetMobileAtoms(atm[icb+1+3*j].UniqueId(),
                                atm[icb+nsatm-1].UniqueId());

}

