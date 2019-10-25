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

#include "Rot.hh"

using namespace UnivConstants;

namespace prf
{
    Rot::Rot() : Update()
    {
        Name("Rot");
        isscu=true;
        isrgu=isbbu=ismcu=false;
        ranges.resize(nranges=1);
        atom_ranges.resize(natomranges=1);
    }

    Rot::~Rot() {}

    void Rot::build_dof_list()
    {
        std::deque<DOF_Info> mydofs;
        for (std::vector<DOF_Info>::iterator i=popl->dof_map().begin();
        i!=popl->dof_map().end();++i) {
            if (i->dof_kind==sidechain_torsion_angle) mydofs.push_back(*i);
        }
        dof_at_site.resize(mydofs.size());
        dof_at_site.assign(mydofs.begin(),mydofs.end());
        site_weight.resize(dof_at_site.size(),1);
    }

    int Rot::perform()
    {
        Update::perform();
        DOF_Info ang=dof_at_site[site];
        thechange[0].before=popl->get_dof(ang.global_index);
        thechange[0].after=twoPi*rnd->shoot();
        // The value "after" could also be chosen from some distribution above
        Ligand * curlg=popl->ligand(ang.group);
        ranges[0].first=ranges[0].second=curlg->UniqueId();

        curlg->rotDof_assign_and_reconstruct(ang.specific_index_in_group,
                                             (thechange[0].after),st_atom,nd_atom);
        atom_ranges[0].first=st_atom;atom_ranges[0].second=nd_atom;
        ++nd_atom;
        thechange[0].info=ang;
        return 0;
    }
}
