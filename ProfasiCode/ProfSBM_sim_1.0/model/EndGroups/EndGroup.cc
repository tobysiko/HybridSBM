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

#include "EndGroup.hh"
using std::deque;
using std::pair;

using namespace prf;
EndGroup::EndGroup(OneLetterCode cod) : Ligand(cod) {}

EndGroup::EndGroup(const EndGroup &geg) : Ligand(geg), bv0(geg.bv0), bv1(geg.bv1),
        rf0(geg.rf0), rf1(geg.rf1), res(geg.res) {}

EndGroup & EndGroup::operator=(const EndGroup &geg)
{
    if (this!=&geg) {
        Ligand::operator=(geg);
        bv0=geg.bv0;
        bv1=geg.bv1;
        rf0=geg.rf0;
        rf1=geg.rf1;
        res=geg.res;
    }

    return *this;
}

EndGroup::~EndGroup() {}

bool EndGroup::pep_bond_link() const
{
    return true;
}

void EndGroup::pep_bond_atoms(int &i1, int &i2)
{
    i1=atom(0).UniqueId();
    i2=atom(1).UniqueId();
}

void EndGroup::Initialize() {}

void EndGroup::Reconstruct() {}

void EndGroup::Write()
{
    prf::cout << "Using EndGroup::Write\n";
}

void EndGroup::WritePDB(int &istatm, char ch_id, int aaindx, FILE * fp)
{
    for (size_t i=0;i<atm.size();++i) {
        WritePDBline(istatm,ch_id,aaindx,i,fp,1);
    }
}

void EndGroup::WritePDB2(int &istatm, char ch_id, int aaindx, FILE * fp)
{
    for (size_t i=0;i<atm.size();++i) {
        if (atm[i].Species() !=hydrogen)
            WritePDBline(istatm,ch_id,aaindx,i,fp,1);
    }

    for (size_t i=0;i<atm.size();++i) {
        if (atm[i].Species() ==hydrogen)
            WritePDBline(istatm,ch_id,aaindx,i,fp,1);
    }
}

int EndGroup::rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1){return 0;}

int EndGroup::rotDof_assign(int il, double mgd) {return 0;}

double EndGroup::get_rotDof(int il) {return 0;}

void EndGroup::LocPairsatRTdof(int i, deque < pair < int, int > >&lcp) {}

void EndGroup::BuildConnections() {}

void EndGroup::Allocate()
{
    Ligand::Allocate();

    for (int i=0;i<mygrp->num_atoms();++i) atm[i].UniqueId(i);
}

void EndGroup::nodes_reconnect(Atom &a1, Atom &a2)
{
    SetRefa(a1,a2);
    if (node.empty()) return;
    int myroot=rf1.UniqueId();
    int n0root=node[0]->atom(1).UniqueId();
    for (size_t i=0;i<node.size();++i) node[i]->AtomOffset(myroot-n0root);
}

Node * EndGroup::node_for_dof(int i)
{
    if (i>=0 and i<nrtdof) return node[i+1];
    //because first node will normally be locked
    else return NULL;
}
