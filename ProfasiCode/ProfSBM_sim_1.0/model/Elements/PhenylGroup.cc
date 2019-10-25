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

#include "PhenylGroup.hh"

using namespace prf;
using std::vector;
using std::deque;
using std::pair;

using namespace prf_utils;

using namespace UnivConstants;
double PhenylGroup::ringbl=1.39;
double PhenylGroup::bXH=1.0;

PhenylGroup::PhenylGroup()
{
    th01=acos(-1.0/3);
    bstm=bprm=1.52;
    nm="Phenyl Group";
    atm.resize(13);   //a reference atm for phi, the "Base",and C6H5
}

PhenylGroup::~PhenylGroup() {}

void PhenylGroup::Initialize()
{
    hexarm.Initialize(bprm,bstm,ringbl,pi-th01,pi/3);
}

void PhenylGroup::Create()
{
    Vector3 v0,v1,v2,org;
    Vector3 hexv0,hexv1,hexv2,rov,hv;
    v0=atm[1].Position()-atm[0].Position();
    hexv0=v1=(org=atm[2].Position())-atm[1].Position();
    v2=v1.cross(v0);
    hexv0.mag(ringbl);
    hexv1=hexarm(v0,v1,v2,phi);
    hexv2=hexarm(v0,v1,v2,phi-pi);
    atm[3].Pos(rov=(org+hexv1));
    atm[4].Pos(rov+hexv0);atm[5].Pos(org+2*hexv0);
    atm[7].Pos(rov=(org+hexv2));atm[6].Pos(rov+hexv0);
    hexv0.mag(bXH);hexv1.mag(bXH);hexv2.mag(bXH);
    atm[8].Pos(atm[3].Position()-hexv2);
    atm[9].Pos(atm[4].Position()+hexv1);
    atm[10].Pos(atm[5].Position()+hexv0);
    atm[11].Pos(atm[6].Position()+hexv2);
    atm[12].Pos(atm[7].Position()-hexv1);
}

void PhenylGroup::ExportConnections(ConnectionsMatrix &aa)
{
    aa.set_connection(atm[0].UniqueId(),atm[1].UniqueId());
    aa.set_connection(atm[0].UniqueId(),atm[2].UniqueId());

    for (unsigned int i=1;i<atm.size();++i)
        for (unsigned int j=i;j<atm.size();++j)
            aa.set_connection(atm[i].UniqueId(),atm[j].UniqueId());

    for (int i=0;i<naltbases;++i) {
        aa.set_connection(altbase[i].UniqueId(),atm[3].UniqueId());
        aa.set_connection(altbase[i].UniqueId(),atm[7].UniqueId());
        aa.set_connection(altbase[i].UniqueId(),atm[5].UniqueId());
        aa.set_connection(altbase[i].UniqueId(),atm[10].UniqueId());
    }
}

void PhenylGroup::LocPairs(deque<pair<int,int> > & lcp)
{
    for (int i=0;i<naltbases;++i) {
        lcp.push_back(std::make_pair(atm[3].UniqueId(),altbase[i].UniqueId()));
        lcp.push_back(std::make_pair(atm[7].UniqueId(),altbase[i].UniqueId()));
    }
}

void PhenylGroup::SetBondLengths(double x0,double x1,double x2,
                                 double x3,double x4)
{
    bprm=x0;bstm=x1;
}

void PhenylGroup::SetBondAngles(double x0,double x1,double x2,double x3)
{
    th01=x0;
}
