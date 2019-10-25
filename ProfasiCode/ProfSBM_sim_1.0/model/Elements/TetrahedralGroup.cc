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

#include "TetrahedralGroup.hh"

using namespace prf;
using std::vector;

using namespace prf_utils;

using namespace UnivConstants;
TetrahedralGroup::TetrahedralGroup()
{
    bstm=bprm=1.52;brnc=1.0;
    th01=th12=acos(-1.0/3);
    nm="Tetrahedral Group";
    atm.resize(6);
}

TetrahedralGroup::TetrahedralGroup(Atom &a0,Atom &a1,Atom &a2,
                                   vector<Atom> &al,int st)
{
    bstm=bprm=1.52;brnc=1.0;
    th01=th12=acos(-1.0/3);
    atm.resize(6);
    AssignAtoms(a0,a1,a2,al,st);
}

TetrahedralGroup::~TetrahedralGroup() {}

void TetrahedralGroup::Initialize() {Initialize(bprm,bstm,brnc,th01,th12);}

void TetrahedralGroup::Initialize(double b0,double b1,double b2,
                                  double th0,double th1)
{
    locate_H.Initialize(b0,b1,b2,pi-th0,pi-th1);
}

void TetrahedralGroup::Create()
{
    Vector3 v0,v1,v2,org;
    v0=atm[1].Position()-atm[0].Position();
    v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);
    atm[3].Pos(org+locate_H(v0,v1,v2,phi));
    atm[4].Pos(org+locate_H(v0,v1,v2,phi+twoPid3));
    atm[5].Pos(org+locate_H(v0,v1,v2,phi-twoPid3));
}

void TetrahedralGroup::SetBondLengths(double x0,double x1,double x2,
                                      double x3,double x4)
{
    bprm=x0;bstm=x1;brnc=x2;
}

void TetrahedralGroup::SetBranchLength(int i, double x)
{
    brnc=x;
}

void TetrahedralGroup::SetBondAngles(double x0,double x1,double x2,double x3)
{
    th01=x0;th12=x1;
}

