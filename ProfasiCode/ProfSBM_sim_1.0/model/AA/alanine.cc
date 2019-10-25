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

#include "alanine.hh"

using namespace UnivConstants;

using namespace prf;
Alanine::Alanine() : AminoAcid(A) {node.resize(1);nnodes=1;}

Alanine::~Alanine() {}

void Alanine::Initialize()
{
    const double bXH=1.0;
    double thH;
    thH=acos(-1.0/3);
    grp.AssignAtoms(Nitrogen(),Calpha(),Cbeta(),atm,ibbca+3);
    grp.Initialize(b[1],bCaCb,bXH,thCb,thH);
    node[0]=&grp;
}
