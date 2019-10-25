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

#include "ATetGroup.hh"

using namespace prf;
using std::vector;

using namespace prf_utils;

using namespace UnivConstants;
ATetGroup::ATetGroup()
{
    bstm=bprm=0;brnca=brncb=brncc=0;
    th01=th1a=th1b=th1c=acos(-1.0/3);
    phab=twoPid3; phac=2*phab;
    atm.resize(6);   //a reference atm for phi, the "Base",and CH3
    nm="ATetGroup";
}

ATetGroup::ATetGroup(Atom &a0,Atom &a1,Atom &a2,vector<Atom> &al,int st)
{
    ATetGroup();
    AssignAtoms(a0,a1,a2,al,st);
    atm.resize(6);   //a reference atm for phi, the "Base",and CH3
}

ATetGroup::~ATetGroup() {}

void ATetGroup::Initialize()
{
    locate_Ha.Initialize(bprm,bstm,brnca,pi-th01,pi-th1a);
    locate_Hb.Initialize(bprm,bstm,brncb,pi-th01,pi-th1b);
    locate_Hc.Initialize(bprm,bstm,brncc,pi-th01,pi-th1c);
}

void ATetGroup::Create()
{
    Vector3 v0,v1,v2,org;
    v0=atm[1].Position()-atm[0].Position();
    v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);
    atm[3].Pos(org+locate_Ha(v0,v1,v2,phi));
    atm[4].Pos(org+locate_Hb(v0,v1,v2,phi+phab));
    atm[5].Pos(org+locate_Hc(v0,v1,v2,phi+phac));
}


void ATetGroup::SetBondLengths(double x0,double x1,double x2,
                               double x3,double x4)
{
    bprm=x0;bstm=x1;brnca=x2;brncb=x3;brncc=x4;
}

void ATetGroup::SetBranchLength(int i, double x)
{
    switch (i%3) {
        case 0: brnca=x; break;
        case 1: brncb=x; break;
        case 2:
        default: brncc=x; break;
    };
}

void ATetGroup::SetBondAngles(double x0,double x1,double x2,double x3)
{
    th01=x0;th1a=x1;th1b=x2;th1c=x3;
}

void ATetGroup::SetRelPhi(double x0,double x1)
{
    phab=x0;phac=x1;
}
