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

#include "HistidineTip.hh"

using namespace prf;
using std::vector;
using std::pair;
using std::deque;

using namespace prf_utils;

using namespace UnivConstants;

HistidineTip::HistidineTip()
{
    th01=113.2*pi/180;th12=126.0*pi/180;
    bprm=bstm=1.52;
    ringbl=1.35;bXH=1.0;
    cu=(cos(0.3*pi)+cos(0.1*pi))/sqrt(2*(1+cos(0.6*pi)));
    cv=0.5/sqrt(2*(1-cos(0.6*pi)));
    cw=0.4*(1+cos(0.3*pi)/(cos(0.3*pi)+cos(0.1*pi)));
    nm="Histidine Tip";
    atm.resize(10);
}

HistidineTip::~HistidineTip() {}

void HistidineTip::Initialize()
{
    double ux;
    pentarm.Initialize(bprm,bstm,ringbl,pi-th01,pi-th12);
    ux=0.4*ringbl*(cos(0.1*pi)+2*cos(0.3*pi));
    cq=(bXH+ux)/ux;
}

void HistidineTip::Create()
{
    Vector3 v0,v1,v2,org;
    Vector3 pv1,pv2,uu,vv;
    v0=atm[1].Position()-atm[0].Position();
    v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);

    pv1=pentarm(v0,v1,v2,phi);
    pv2=pentarm(v0,v1,v2,phi-pi);
    atm[3].Pos(org+pv1);
    atm[6].Pos(org+pv2);
    uu=cu*(pv1+pv2);vv=cv*(pv1-pv2);
    atm[4].Pos(org+uu+vv);atm[5].Pos(org+uu-vv);
    org+=cw*uu;
    atm[7].Pos(org+cq*(atm[3].Position()-org));
    atm[8].Pos(org+cq*(atm[5].Position()-org));
    atm[9].Pos(org+cq*(atm[6].Position()-org));
}

void HistidineTip::LocPairs(deque<pair<int,int> > & lcp)
{
    for (int i=0;i<naltbases;++i) {
        lcp.push_back(std::make_pair(atm[3].UniqueId(),altbase[i].UniqueId()));
        lcp.push_back(std::make_pair(atm[6].UniqueId(),altbase[i].UniqueId()));
    }
}

void HistidineTip::ExportConnections(ConnectionsMatrix &aa)
{
    aa.set_connection(atm[0].UniqueId(),atm[1].UniqueId());
    aa.set_connection(atm[0].UniqueId(),atm[2].UniqueId());

    for (unsigned int i=1;i<atm.size();++i)
        for (unsigned int j=i;j<atm.size();++j)
            aa.set_connection(atm[i].UniqueId(),atm[j].UniqueId());

    for (int i=0;i<naltbases;++i) {
        aa.set_connection(altbase[i].UniqueId(),atm[3].UniqueId());
        aa.set_connection(altbase[i].UniqueId(),atm[6].UniqueId());
    }
}

void HistidineTip::SetBondLengths(double x0,double x1,double x2,
                                  double x3,double x4)
{
    bprm=x0;bstm=x1;ringbl=x2;
}

void HistidineTip::SetBondLengths(double x0,double x1,double x2)
{
    bprm=x0;bstm=x1;ringbl=x2;
}
