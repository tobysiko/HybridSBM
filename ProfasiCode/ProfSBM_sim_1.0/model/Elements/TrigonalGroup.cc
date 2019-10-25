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

#include "TrigonalGroup.hh"

using namespace prf;

using namespace prf_utils;

using namespace UnivConstants;

TrigonalGroup::TrigonalGroup()
{
    bstm=bprm=1.52;brnc=1.0;
    th0=acos(-1.0/3);th1=twoPid3;
    atm.clear();
    atm.resize(5);   //a reference atm for phi, the "Base",and NH2
    nm="Trigonal Group";
}

TrigonalGroup::~TrigonalGroup() {}

void TrigonalGroup::Initialize()
{
    locate_H.Initialize(bprm,bstm,brnc,pi-th0,pi-th1);
}

void TrigonalGroup::Create()
{
    Vector3 v0,v1,v2,org;
    v0=atm[1].Position()-atm[0].Position();
    v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);
    atm[3].Pos(org+locate_H(v0,v1,v2,phi));
    atm[4].Pos(org+locate_H(v0,v1,v2,phi+pi));
}

void TrigonalGroup::SetBondLengths(double x0,double x1,double x2,
                                   double x3,double x4)
{
    bprm=x0;bstm=x1;brnc=x2;
}

void TrigonalGroup::SetBranchLength(int i, double x)
{
    brnc=x;
}

void TrigonalGroup::SetBondAngles(double x0,double x1,double x2,double x3)
{
    th0=x0;th1=x1;
}
