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

#include "JacobianCarTor.hh"

JacobianCarTor::JacobianCarTor() : popl(NULL), nat(0) {}

JacobianCarTor::~JacobianCarTor() {}

void JacobianCarTor::set_population(prf::Population *popu)
{
    if (popu==NULL) return;
    popl=popu;
    nat=popl->NumberOfAtoms();
    Atom a1,a2;
    int nalldof=popl->n_dof();
    axis.resize(nalldof,std::make_pair(-2,-2));
    chlist.resize(nalldof,std::make_pair(0,-1));
    chstrt.resize(popl->NumberOfChains(),0);

    for (size_t idof=0;idof<axis.size();++idof) {
        DOF_Info dof=popl->get_dof_info(idof);
        if (dof.dof_kind==backbone_torsion_angle) {
            int idir=popl->Chain(dof.chain)
                     ->BBAxisAtoms(a1,a2,dof.specific_index_in_chain);
            axis[idof].first=a1.UniqueId();
            axis[idof].second=a2.UniqueId();
            chlist[idof].first=a2.UniqueId();
            if (idir==0) {
                // Origin for the rotation is at an atom in the second half of
                // the chain.
                chlist[idof].second=popl->Chain(dof.chain)->last_atom();
            } else {
                // Origin for the rotation is in the first half of the chain.
                // Notice that the range goes from a higher to a lower value.
                // This encodes the instruction that the contributions to the
                // Jacobian must be multiplied by -1.
                chlist[idof].second=popl->Chain(dof.chain)->first_atom();
            }
        } else if (dof.dof_kind==sidechain_torsion_angle) {
            prf::Ligand *lg=popl->ligand(dof.group);
            lg->get_rotDofAxis(dof.specific_index_in_group,a1,a2);
            axis[idof].first=a1.UniqueId();
            axis[idof].second=a2.UniqueId();
            int i1,i2;
            lg->node_for_dof(dof.specific_index_in_group)->MobileAtoms(i1,i2);
            chlist[idof].first=i1;
            chlist[idof].second=i2;
        } else if (dof.dof_kind==rigid_body_xyz) {
            axis[idof].first=-1;
            axis[idof].second=dof.chain;
            if (dof.index_in_chain==0) chstrt[dof.chain]=idof;
            chlist[idof].first=popl->Chain(dof.chain)->first_atom();
            chlist[idof].second=popl->Chain(dof.chain)->last_atom();
        }
    }
    allocate(popl->n_dof(),3*nat);
}

void JacobianCarTor::set_dofmap(std::vector<int> &indxs)
{
    indexmap=indxs;
    allocate(indexmap.size(),3*nat);
}

void JacobianCarTor::refresh()
{
    prf::Vector3 org,a,r,dr;
    v=0;
    for (size_t i=0;i<indexmap.size();++i) {
        std::pair<int,int> ax=axis[indexmap[i]];
        if (ax.first>=0 and ax.first<nat) {
            org=popl->atom(ax.first).Position();
            a=popl->atom(ax.second).Position()-org;
            a.mag(1);
            std::pair<int,int> ch=chlist[indexmap[i]];
            if (ch.second>ch.first) {
                for (int j=ch.first;j<=ch.second;++j) {
                    r=popl->atom(j).Position();
                    dr=a.cross(r-org);
                    for (int k=0;k<3;++k) {
                        set(i,3*j+k,dr[k]);
                    }
                }
            } else {
                for (int j=ch.second;j<=ch.first;++j) {
                    r=popl->atom(j).Position();
                    dr=(r-org).cross(a);
                    for (int k=0;k<3;++k) {
                        set(i,3*j+k,dr[k]);
                    }
                }
            }
        } else if (ax.first==-1) {
            int ofst=(indexmap[i]-chstrt[ax.second])%3;
            for (int j=popl->Chain(ax.second)->begin_atom();
            j<popl->Chain(ax.second)->end_atom();++j) {
                set(i,3*j+ofst,1);
            }
        }
    }
}
