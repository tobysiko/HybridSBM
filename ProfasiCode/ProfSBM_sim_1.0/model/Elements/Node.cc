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

#include "Node.hh"
#include "../Aux/Constants.hh"
using std::vector;
using std::deque;
using std::pair;
using std::make_pair;

using namespace prf;

Node::Node() : phi(UnivConstants::pi),rotlocked(0),upst(0),upnd(0)
{
    subnode[0]=subnode[1]=NULL;nsubnd=0;naltbases=0;Name("Node Base");
}

Node::~Node() {}

Node::Node(const Node &gnd) : Named(gnd), phi(gnd.phi), rotlocked(gnd.rotlocked),
        atm(gnd.atm), naltbases(gnd.naltbases), nsubnd(gnd.nsubnd), upst(gnd.upst),
        upnd(gnd.upnd)
{
    for (int i=0;i<3;++i) altbase[i]=gnd.altbase[i];

    subnode[0]=gnd.subnode[0]; subnode[1]=gnd.subnode[1];
}

Node &Node::operator=(const Node &gnd)
{
    if (this!=&gnd) {
        phi=gnd.phi;
        rotlocked=gnd.rotlocked;
        atm=gnd.atm;
        naltbases=gnd.naltbases;
        nsubnd=gnd.nsubnd;
        upst=gnd.upst;
        upnd=gnd.upnd;

        for (int i=0;i<3;++i) altbase[i]=gnd.altbase[i];

        subnode[0]=gnd.subnode[0]; subnode[1]=gnd.subnode[1];
    }

    return *this;
}

void Node::Initialize() {}

void Node::AssignAtoms(Atom &a0, Atom &a1,Atom &a2,vector<Atom> &att,int st)
{
    atm[0]=a0;atm[1]=a1;atm[2]=a2;

    for (size_t i=3;i<atm.size();++i) atm[i]=att[st+i-3];

    SetMobileAtoms(atm[3].UniqueId(),atm.back().UniqueId());

}

void Node::AssignAtoms(Atom &a0, Atom &a1, Atom &a2, Atom &a3)
{
    atm[0]=a0;atm[1]=a1;atm[2]=a2;

    if (atm.size()>3) atm[3]=a3;

    SetMobileAtoms(atm[3].UniqueId(),atm.back().UniqueId());
}

void Node::AssignAtoms(Atom &a0, Atom &a1,Atom &a2,
                       Atom &a3,Atom &a4)
{
    atm[0]=a0;atm[1]=a1;atm[2]=a2;

    if (atm.size()>3) atm[3]=a3;

    if (atm.size()>4) atm[4]=a4;

    SetMobileAtoms(atm[3].UniqueId(),atm.back().UniqueId());
}

void Node::AssignAtoms(Atom &a0, Atom &a1,Atom &a2,
                       Atom &a3,Atom &a4,Atom &a5)
{
    atm[0]=a0;atm[1]=a1;atm[2]=a2;

    if (atm.size()>3) atm[3]=a3;

    if (atm.size()>4) atm[4]=a4;

    if (atm.size()>5) atm[5]=a5;

    SetMobileAtoms(atm[3].UniqueId(),atm.back().UniqueId());
}

void Node::Create()
{
    cerr<<"Create called for base node class\n";
}

void Node::ReCreate()
{
    Create();

    for (int i=0;i<nsubnd;++i) {
        subnode[i]->ReCreate();
    }
}

void Node::AssignPhi(double x)
{
    phi=x;ReCreate();
}

void Node::RevertPhi(double x) {phi=x;}

void Node::AddSubnode(Node *nd) {if (nd!=NULL) {subnode[nsubnd++]=nd;}}

void Node::SetBase(Atom &a0,Atom &a1, Atom &a2)
{
    atm[0]=a0;atm[1]=a1;atm[2]=a2;
}

void Node::MobileAtoms(int &i, int &j) {i=upst;j=upnd;}

void Node::SetMobileAtoms(int i, int j)
{
    if (j>=i) {
        upst=i;upnd=j;
    } else {
        upst=j;upnd=i;
    }
}

void Node::AtomOffset(int i)
{
    for (unsigned int j=0;j<atm.size(); ++j) {
        atm[j].UniqueId(atm[j].UniqueId()+i);
    }

    for (int j=0;j<naltbases;++j) {
        altbase[j].UniqueId(altbase[j].UniqueId()+i);
    }

    upst+=i;upnd+=i;
}

void Node::ExportConnections(ConnectionsMatrix &aa)
{
    aa.set_connection(atm[0].UniqueId(),atm[1].UniqueId());
    aa.set_connection(atm[0].UniqueId(),atm[2].UniqueId());

    for (unsigned int i=1;i<atm.size();++i)
        for (unsigned int j=i;j<atm.size();++j)
            aa.set_connection(atm[i].UniqueId(),atm[j].UniqueId());

    for (int si=0;si<nsubnd;++si) subnode[si]->ExportConnections(aa);

    for (int i=0;i<naltbases;++i) {
        for (unsigned int j=0;j<atm.size();++j) {
            aa.set_connection(altbase[i].UniqueId(),atm[j].UniqueId());
        }
    }
}

void Node::LocPairs(deque<pair<int,int> > & lcp)
{
    lcp.clear();

    if (rotlocked!=1) {
        for (size_t i=3;i<atm.size();++i) {
            for (int j=0;j<naltbases;++j) {
                lcp.push_back(make_pair(atm[i].UniqueId(),
                                        altbase[j].UniqueId()));
            }
        }
    }
}

void Node::BuildConnections()
{
    for (int si=0;si<nsubnd;++si) {
        subnode[si]->SetBase(atm[1]);

        for (unsigned int ab=3;ab<atm.size();++ab) {
            if (atm[ab].UniqueId()!=(subnode[si]->ATOM(2).UniqueId()))
                subnode[si]->SetBase(atm[ab]);
        }

        subnode[si]->BuildConnections();
    }
}

void Node::SetBase(Atom & ap)
{
    if (naltbases<=2) altbase[naltbases++]=ap; else altbase[2]=ap;
}

void Node::set_bases(Atom &b1, Atom &b2, Atom &b3)
{
    altbase[0]=b1;
    altbase[1]=b2;
    altbase[2]=b3;
    naltbases=3;
}

void Node::set_bases(Atom &b1, Atom &b2)
{
    altbase[0]=b1;
    altbase[1]=b2;
    naltbases=2;
}

void Node::SetBondLengths(double x0,double x1,double x2,double x3,double x4){}

void Node::SetBondAngles(double x0,double x1,double x2,double x3) {}

void Node::SetRelPhi(double x0,double x1){}

void Node::SetBranchLength(int i, double x) {}
