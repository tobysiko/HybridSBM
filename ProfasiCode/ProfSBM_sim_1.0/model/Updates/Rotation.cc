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

#include "Rotation.hh"

using namespace UnivConstants;

namespace prf
{
    Rotation::Rotation() : Update()
    {
        Name("Rotation");scl=0.1*pi;
        isrgu=true;
        isscu=isbbu=ismcu=false;
        ranges.resize(1);
        atom_ranges.resize(1);
        nranges=1;natomranges=1;
        nchanges=9;
        thechange.resize(9);
    }

    Rotation::~Rotation() {}

    void Rotation::setScale(double gscl)
    {
        scl=gscl;
        Logger(5)<<Name()<<"> setting angular step size to "<<scl<<"\n";
    }

    void Rotation::build_dof_list()
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

    int Rotation::perform()
    {
        double theta,phi;
        Update::perform();
        int chid=dof_at_site[site].chain;
        int sitebase=site-site%9;
        theta=acos(2*rnd->shoot()-1); // distribute theta like sin(theta)
        phi=twoPi*rnd->shoot();         // to avoid over density at poles
        Vector3 rax(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        double rotangl=scl*rnd->shoot();
        atom_ranges[0].first=st_atom=popl->Chain(chid)->begin_atom();
        atom_ranges[0].second=(nd_atom=popl->Chain(chid)->end_atom())-1;
        ranges[0].first=popl->Chain(chid)->first_ligand()->UniqueId();
        ranges[0].second=popl->Chain(chid)->last_ligand()->UniqueId();
        AtomCoordinates midatomloc((st_atom+nd_atom)/2);
        Vector3 org(midatomloc.value());
        //cout <<"Rotation of atoms from "<<st_atom<<" to "<<nd_atom<<"\n";
        //cout <<"by ammount "<<rotangl<<"\n";
        for (size_t i=0;i<9;++i) {
            thechange[i].info=dof_at_site[sitebase+i];
            thechange[i].before=popl->get_dof(dof_at_site[sitebase+i]);
        }

        AtomCoordinates::BlockRotate(rotangl,st_atom,nd_atom,org,rax);

        for (size_t i=0;i<9;++i) {
            thechange[i].after=popl->get_dof(dof_at_site[sitebase+i]);
        }
        return 0;
    }

}
