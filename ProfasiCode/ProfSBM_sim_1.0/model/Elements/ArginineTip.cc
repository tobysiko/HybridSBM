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

#include "ArginineTip.hh"

using namespace prf;
using std::vector;
using std::deque;
using std::pair;

using namespace prf_utils;

using namespace UnivConstants;
ArginineTip::ArginineTip()
{
    bprm=1.52;bstm=1.46;brncnc=1.33;brnccn=1.33;brncnh=1.0;
    th01=111.5*pi/180;th12=124.4*pi/180;th23=th34=120.0*pi/180;
    nhdcn=brncnh/brnccn;
    nhdnc=brncnh/brncnc;
    nm="ArginineTip";
    atm.resize(11);
}

ArginineTip::ArginineTip(Atom &a0,Atom &a1,Atom &a2,vector<Atom> &al,int st)
{
    ArginineTip();
    AssignAtoms(a0,a1,a2,al,st);
}

ArginineTip::~ArginineTip() {}

void ArginineTip::Initialize()
{
    locate_C.Initialize(bprm,bstm,brncnc,pi-th01,pi-th12);
    locate_N.Initialize(bstm,brncnc,brnccn,pi-th12,pi-th23);
}

void ArginineTip::Create()
{
    Vector3 v0,v1,v2,v3,org,orgn,orgc,vc,vn1,vn2;
    v0=atm[1].Position()-atm[0].Position();
    v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);
    vc=locate_C(v0,v1,v2,phi);
    v3=vc.cross(v1);
    vn1=locate_N(v1,vc,v3,pi);vn2=locate_N(v1,vc,v3,0);
    atm[4].Pos(orgc=org+vc);
    atm[5].Pos(orgn=(orgc+vn1));
    atm[6].Pos(orgn+nhdnc*vc);
    atm[7].Pos(orgn-nhdcn*vn2);
    atm[8].Pos(orgn=orgc+vn2);
    atm[9].Pos(orgn-nhdcn*vn1);
    atm[10].Pos(orgn+nhdnc*vc);
    atm[3].Pos(org-nhdcn*vn2);
}

void ArginineTip::LocPairs(deque<pair<int,int> > & lcp)
{
    for (int i=0;i<naltbases;++i) {
        lcp.push_back(std::make_pair(atm[3].UniqueId(),altbase[i].UniqueId()));
        lcp.push_back(std::make_pair(atm[4].UniqueId(),altbase[i].UniqueId()));
    }
}

void ArginineTip::ExportConnections(ConnectionsMatrix &aa)
{
    aa.set_connection(atm[0].UniqueId(),atm[1].UniqueId());
    aa.set_connection(atm[0].UniqueId(),atm[2].UniqueId());

    for (unsigned int i=1;i<atm.size();++i)
        for (unsigned int j=i;j<atm.size();++j)
            aa.set_connection(atm[i].UniqueId(),atm[j].UniqueId());

    for (int i=0;i<naltbases;++i) {
        aa.set_connection(altbase[i].UniqueId(),atm[3].UniqueId());
        aa.set_connection(altbase[i].UniqueId(),atm[4].UniqueId());
    }
}
