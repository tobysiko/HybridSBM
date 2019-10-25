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

#include "DOF_Index.hh"
#include "DOF_Info.hh"

using namespace prf;
key_pair::key_pair() : k1(-1), k2(-1) {}
key_pair::key_pair(int i1,int i2) : k1(i1), k2(i2) {}
key_pair::~key_pair() {}
key_pair::key_pair(const key_pair &k) : k1(k.k1), k2(k.k2) {}
key_pair &key_pair::operator=(const key_pair &k)
                             {
    if (this!=&k) {
        k1=k.k1;
        k2=k.k2;
    }
    return *this;
}

bool key_pair::operator==(const key_pair &k) const
{
    return (k1==k.k1 and k2==k.k2);
}

bool key_pair::operator<(const key_pair &k) const
{
    return k1<k.k1 or (k1==k.k1 && k2<k.k2);
}

key_tripplet::key_tripplet() : k1(-1),k2(-1), k3(-1) {}
key_tripplet::key_tripplet(int i1,int i2, int i3) : k1(i1), k2(i2), k3(i3) {}
key_tripplet::~key_tripplet() {}
key_tripplet::key_tripplet(const key_tripplet &k) : k1(k.k1),k2(k.k2),k3(k.k3)
{}
key_tripplet &key_tripplet::operator=(const key_tripplet &k)
                                     {
    if (this!=&k) {
        k1=k.k1,
        k2=k.k2;
        k3=k.k3;
    }
    return *this;
}
bool key_tripplet::operator==(const key_tripplet &k) const
{
    return (k1==k.k1 and k2==k.k2 and k3==k.k3);
}

bool key_tripplet::operator<(const key_tripplet &k) const
{
    return (k1<k.k1 or (k1==k.k1 and (k2<k.k2 or (k2==k.k2 and k3<k.k3))));
}


DOF_Index::DOF_Index() {}

DOF_Index::~DOF_Index() {}

int DOF_Index::init(std::vector<DOF_Info> &dofs)
{
    for (size_t i=0;i<dofs.size();++i) {
        int uid=dofs[i].global_index;
        key_pair t2(dofs[i].chain, dofs[i].index_in_chain);
        ch_chind[t2]=uid;
        t2.assign(dofs[i].dof_kind, dofs[i].specific_global_index);
        tp_tpind[t2]=uid;

        key_tripplet t3(dofs[i].chain,dofs[i].dof_kind,dofs[i].specific_index_in_chain);
        ch_t_chtind[t3]=uid;
        t3.assign(dofs[i].dof_kind,dofs[i].group,dofs[i].specific_index_in_group);
        rs_t_rstind[t3]=uid;
    }
    return 1;
}

void DOF_Index::reset()
{
    ch_chind.clear();
    tp_tpind.clear();
    ch_t_chtind.clear();
    rs_t_rstind.clear();
}

int DOF_Index::get_uid_from_chain_index(int ch, int chind)
{
    key_pair kk(ch,chind);
    std::map<key_pair,int>::iterator it=ch_chind.find(kk);
    if (it!=ch_chind.end()) return it->second;
    else return -1;
}

int DOF_Index::get_uid_from_type_index(int tp, int tpind)
{
    key_pair kk(tp,tpind);
    std::map<key_pair,int>::iterator it=tp_tpind.find(kk);
    if (it!=tp_tpind.end()) return it->second;
    else return -1;
}

int DOF_Index::get_uid_from_chain_type(int ch, int tp, int chtpind)
{
    key_tripplet kkk(ch,tp,chtpind);
    std::map<key_tripplet,int>::iterator it=ch_t_chtind.find(kkk);
    if (it!=ch_t_chtind.end()) return it->second;
    else return -1;
}

int DOF_Index::get_uid_from_residue_type(int rs, int tp, int rstpind)
{
    key_tripplet kkk(tp,rs,rstpind);
    std::map<key_tripplet,int>::iterator it=rs_t_rstind.find(kkk);
    if (it!=rs_t_rstind.end()) return it->second;
    else return -1;
}
