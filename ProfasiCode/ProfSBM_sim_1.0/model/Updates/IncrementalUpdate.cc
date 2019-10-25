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

#include "IncrementalUpdate.hh"

IncrementalUpdate::IncrementalUpdate()
{
    isscu=false;
    isrgu=isbbu=ismcu=islocal=issequencial=false;
    st_fl=nd_fl=st_atom=nd_atom=0;
    state=0;
    itmp=0;
    nchanges=1;
    dirn=0;
    thechange.resize(1);
    ranges.resize(2);
    atom_ranges.resize(1);
    nranges=natomranges=0;
    Name("IncrementalUpdate");
    rnd=NULL;
}

IncrementalUpdate::~IncrementalUpdate() {}

void IncrementalUpdate::pick_dof(int i)
{
    dof=popl->get_dof_info(i);
    if (dof.global_index<0) return;
    if (dof.dof_kind==backbone_torsion_angle) {
        isscu=ismcu=isrgu=false;
        isbbu=true;
        nranges=2;
        natomranges=1;
    } else if (dof.dof_kind==sidechain_torsion_angle) {
        isbbu=ismcu=isrgu=false;
        isscu=true;
        nranges=1;
        natomranges=1;
    }
}

int IncrementalUpdate::perform()
{
    state=1;
    thechange[0].before=popl->get_dof(dof.global_index);
    thechange[0].after=thechange[0].before+scl;
    if (isscu) {
        Ligand * curlg=popl->ligand(dof.group);
        ranges[0].first=
                ranges[0].second=
                curlg->UniqueId()-popl->chain_start(dof.chain);
        curlg->rotDof_assign_and_reconstruct(dof.specific_index_in_group,
                                             (thechange[0].after),st_atom,nd_atom);
        atom_ranges[0].first=st_atom;atom_ranges[0].second=nd_atom;
        ++nd_atom;
        thechange[0].info=dof;
    } else if (isbbu) {
        Atom ax0, ax1;
        int chid=dof.chain;
        dirn = popl->Chain(chid)->BBAxisAtoms(ax0, ax1,
                                              dof.specific_index_in_chain);
        int upsch=dof.group;
        int upaano = upsch-popl->Chain(chid)->first_AA()->UniqueId();

        if (dirn == 0) {
            ranges[0].first=ranges[0].second=upsch;
            ranges[1].first=upsch+1;
            ranges[1].second=popl->Chain(chid)->last_ligand()->UniqueId();
            if (ranges[1].second<ranges[1].first) nranges=1; else nranges=2;

            atom_ranges[0].first=st_atom = ax1.UniqueId();
            atom_ranges[0].second=(nd_atom= popl->Chain(chid)->end_atom())-1;
        } else {
            ranges[0].first=popl->Chain(chid)->first_ligand()->UniqueId();
            ranges[0].second=upsch-1;
            if (ranges[0].second<ranges[0].first) nranges=1; else nranges=2;

            ranges[nranges-1].first=ranges[nranges-1].second=upsch;
            atom_ranges[0].first=st_atom = popl->Chain(chid)->begin_atom();
            atom_ranges[0].second=(nd_atom = ax1.UniqueId())-1;
        }

        if (dirn == 0) {
            AtomCoordinates::BlockRotate(scl, st_atom, nd_atom,
                                         ax0.Pos(), ax1.Pos());
            popl->Chain(chid)->backbone()->freconst_bond_vectors(3*upaano);
        } else {
            AtomCoordinates::BlockRotate(scl, st_atom, nd_atom,
                                         ax1.Pos(), ax0.Pos());
            popl->Chain(chid)->backbone()->breconst_bond_vectors(3*upaano+2);
        }

        popl->Chain(chid)->BBdof(dof.specific_index_in_chain,
                                 thechange[0].after);
        thechange[0].info=dof;
    }
    return 0;
}

int IncrementalUpdate::revert()
{
    if (Update::revert()) {
        if (isbbu) {
            int upch=thechange[0].info.chain;
            int upsch=thechange[0].info.group-popl->Chain(upch)->first_AA()->UniqueId();

            if (dirn == 0) {
                popl->Chain(upch)->backbone()->freconst_bond_vectors(3*upsch);
            } else {
                popl->Chain(upch)->backbone()->breconst_bond_vectors(3*upsch+2);
            }
        }
        return 1;
    } else
        return 0;
}
