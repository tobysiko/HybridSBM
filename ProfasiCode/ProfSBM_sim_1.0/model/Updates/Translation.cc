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

#include "Translation.hh"

using namespace UnivConstants;

namespace prf
{
    Translation::Translation() : Update()
    {
        Name("Translation");
        scl=0.1*AtomCoordinates::boxL();
        isrgu=true;
        isscu=isbbu=ismcu=false;
        ranges.resize(1);
        atom_ranges.resize(1);
        nranges=1;
        natomranges=1;
        nchanges=9;
        thechange.resize(9);
    }

    Translation::~Translation() {}

    void Translation::setScale(double gscl)
    {
        scl=gscl;
        Logger(5)<<Name()<<"> setting step size to "<<scl<<"\n";
    }

    void Translation::build_dof_list()
    {
        std::deque<DOF_Info> mydofs;
        for (std::vector<DOF_Info>::iterator i=popl->dof_map().begin();
        i!=popl->dof_map().end();++i) {
            if (i->dof_kind==rigid_body_xyz) mydofs.push_back(*i);
        }
        dof_at_site.resize(mydofs.size());
        dof_at_site.assign(mydofs.begin(),mydofs.end());
        site_weight.resize(dof_at_site.size(),1);
    }

    int Translation::perform()
    {
        double theta,phi;
        Update::perform();
        int chid=dof_at_site[site].chain;
        int sitebase=site-site%9;
        theta=pi*rnd->shoot();
        phi=twoPi*rnd->shoot();
        Vector3 trv(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        trv*=(rnd->shoot()*scl);
        atom_ranges[0].first=st_atom=popl->Chain(chid)->begin_atom();
        atom_ranges[0].second=(nd_atom=popl->Chain(chid)->end_atom())-1;
        ranges[0].first=popl->Chain(chid)->first_ligand()->UniqueId();
        ranges[0].second=popl->Chain(chid)->last_ligand()->UniqueId();
        //cout <<"Translation of atoms from "<<st_atom<<" to "<<nd_atom<<"\n";
        for (size_t i=0;i<9;++i) {
            thechange[i].info=dof_at_site[sitebase+i];
            thechange[i].before=popl->get_dof(dof_at_site[sitebase+i]);
            thechange[i].after=thechange[i].before+trv[i%3];
        }
        AtomCoordinates::BlockTranslate(trv,st_atom,nd_atom);
        return 0;
    }

}
