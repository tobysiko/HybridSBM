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

#include "ATriGroup.hh"

using namespace prf;

using namespace prf_utils;

using namespace UnivConstants;

ATriGroup::ATriGroup()
{
    bstm=bprm=1.52;brnca=brncb=1.0;phba=pi;
    th0=acos(-1.0/3);th1a=th1b=twoPid3;
    atm.resize(5);   //a reference atm for phi, the "Base",and NH2
    nm="ATriGroup";
}

ATriGroup::~ATriGroup() {}

void ATriGroup::Initialize()
{
    locate_Ha.Initialize(bprm,bstm,brnca,pi-th0,pi-th1a);
    locate_Hb.Initialize(bprm,bstm,brncb,pi-th0,pi-th1b);
}

void ATriGroup::Create()
{
    Vector3 v0,v1,v2,org;
    v0=atm[1].Position()-atm[0].Position();
    v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);
    atm[3].Pos(org+locate_Ha(v0,v1,v2,phi));
    atm[4].Pos(org+locate_Hb(v0,v1,v2,phi+phba));
}

void ATriGroup::SetBondLengths(double x0,double x1,double x2,
                               double x3,double x4)
{
    bprm=x0;bstm=x1;brnca=x2;brncb=x3;
}

void ATriGroup::SetBranchLength(int i, double x)
{
    if (i%2==0) brnca=x; else brncb=x;
}

void ATriGroup::SetBondAngles(double x0,double x1,double x2,double x3)
{
    th0=x0;th1a=x1;th1b=x2;
}

void ATriGroup::SetRelPhi(double x0,double x1)
{
    phba=x0;
}
