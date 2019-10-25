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

#include "Pivot.hh"

using namespace UnivConstants;

namespace prf
{
    Pivot::Pivot():Update()
    {
        Name("Pivot");
        isbbu=true;
        isscu=isrgu=ismcu=false;
        ranges.resize(2);
        atom_ranges.resize(1);
        nranges=2;natomranges=1;
        dirn = 0;
        scl=twoPi;
    }

    Pivot::~Pivot() {}

    void Pivot::build_dof_list()
    {
        std::deque<DOF_Info> mydofs;
        for (std::vector<DOF_Info>::iterator i=popl->dof_map().begin();
        i!=popl->dof_map().end();++i) {
            if (i->dof_kind==backbone_torsion_angle) mydofs.push_back(*i);
        }
        dof_at_site.resize(mydofs.size());
        dof_at_site.assign(mydofs.begin(),mydofs.end());
        site_weight.resize(dof_at_site.size(),1);
    }

    int Pivot::perform()
    {
        Update::perform();
        DOF_Info ang=dof_at_site[site];
        thechange[0].before=popl->get_dof(ang.global_index);
        double rotangl=scl*(rnd->shoot()-0.5);
        thechange[0].after=thechange[0].before+rotangl;

        Atom ax0, ax1;
        int chid=ang.chain;
        dirn = popl->Chain(chid)->BBAxisAtoms(ax0, ax1,
                                              ang.specific_index_in_chain);
        int upsch=ang.group;
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
            AtomCoordinates::BlockRotate(rotangl, st_atom, nd_atom,
                                         ax0.Pos(), ax1.Pos());
            popl->Chain(chid)->backbone()->freconst_bond_vectors(3*upaano);
        } else {
            AtomCoordinates::BlockRotate(rotangl, st_atom, nd_atom,
                                         ax1.Pos(), ax0.Pos());
            popl->Chain(chid)->backbone()->breconst_bond_vectors(3*upaano+2);
        }

        popl->Chain(chid)->BBdof(ang.specific_index_in_chain,
                                 thechange[0].after);
        thechange[0].info=ang;
        return 0;
    }

    int Pivot::revert()
    {
        if (Update::revert()) {
            int upch=thechange[0].info.chain;
            int upsch=thechange[0].info.group-popl->Chain(upch)->first_AA()->UniqueId();

            if (dirn == 0) {
                popl->Chain(upch)->backbone()->freconst_bond_vectors(3*upsch);
            } else {
                popl->Chain(upch)->backbone()->breconst_bond_vectors(3*upsch+2);
            }

            return 1;
        } else
            return 0;
    }
}

