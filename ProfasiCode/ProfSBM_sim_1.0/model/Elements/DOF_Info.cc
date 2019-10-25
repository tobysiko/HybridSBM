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

#include "DOF_Info.hh"
#include <cstdio>

using namespace prf;

DOF_Info::DOF_Info() : dof_kind(rigid_body_xyz), chain(0), group(0),
global_index(0), index_in_chain(0), specific_index_in_group(0),
specific_global_index(0), specific_index_in_chain(0) {}

DOF_Info::~DOF_Info() {}

DOF_Info::DOF_Info(const DOF_Info &d) : dof_kind(d.dof_kind), chain(d.chain),
group(d.group),global_index(d.global_index),
index_in_chain(d.index_in_chain),
specific_index_in_group(d.specific_index_in_group),
specific_global_index(d.specific_global_index),
specific_index_in_chain(d.specific_index_in_chain) {}

DOF_Info & DOF_Info::operator=(const DOF_Info &d)
{
    if (this!=&d) {
        dof_kind=d.dof_kind;
        chain=d.chain;
        group=d.group;
        global_index=d.global_index;
        index_in_chain=d.index_in_chain;
        specific_index_in_group=d.specific_index_in_group;
        specific_global_index=d.specific_global_index;
        specific_index_in_chain=d.specific_index_in_chain;
    }
    return *this;
}

void DOF_Info::assign(int kind, int ch_no, int grp_no, int glob_idx,
                      int idx_in_chn, int sp_idx_in_grp, int sp_glob_idx,
                      int sp_idx_in_chn)
{
    dof_kind=kind;
    chain=ch_no;
    group=grp_no;
    global_index=glob_idx;
    index_in_chain=idx_in_chn;
    specific_index_in_group=sp_idx_in_grp;
    specific_global_index=sp_glob_idx;
    specific_index_in_chain=sp_idx_in_chn;
}

std::string DOF_Info::str() const
{
    char nm[64];
    char tp='r';
    if (dof_kind==backbone_torsion_angle) tp='b';
    else if (dof_kind==sidechain_torsion_angle) tp='s';
    sprintf(nm,"%d::%c:%d",chain,tp,specific_index_in_chain);
    return std::string(nm);
}
