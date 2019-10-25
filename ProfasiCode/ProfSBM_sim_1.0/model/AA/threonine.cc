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

#include "threonine.hh"

using namespace UnivConstants;

using namespace prf;
Threonine::Threonine() :AminoAcid(T)
{
    node.resize(nnodes = 3);
}

Threonine::~Threonine()
{
}

void Threonine::Initialize()
{
    const double c1 = 1.43, c2 = 1.52, bXH = 1.0;
    double th1, th2, thH, thHg1;
    th1 = 108.6 * pi / 180;
    th2 = 111.5 * pi / 180;
    thH = 109.3 * pi / 180;
    thHg1 = 108.0 * pi / 180;
    g0.AssignAtoms(Nitrogen(), Calpha(), Cbeta(),
                   at(" OG1"), at(" HB "), at(" CG2"));
    //The oxygen comes first in the conventional ordering as in the PDB files,
    //which is not the chiral order. Thus the peculiar placement of the HB above.
    g1.AssignAtoms(Calpha(), Cbeta(), at(" OG1"), at(" HG1"));
    g2.AssignAtoms(Calpha(), Cbeta(), at(" CG2"), atm, mygrp->index("1HG2"));
    g0.SetTheta(thCb, th1, thH, th2);
    g1.SetTheta(th1, thHg1);
    g2.SetRootAngle(th2);
    g0.SetBondLengths(b[1], bCaCb, c1, bXH, c2);
    g1.SetBondLengths(bCaCb, c1, bXH);
    g2.SetBondLengths(bCaCb, c2, bXH);
    node[0] = &g0;
    node[1] = &g1;
    node[2] = &g2;
    node[0]->AddSubnode(node[1]);
    node[0]->AddSubnode(node[2]);

    for (int j = 0; j < 3; ++j)
        node[j]->Initialize();

    node[0]->SetMobileAtoms(at(" HB ").UniqueId(),
                            at(" HG1").UniqueId());
}

