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

#include "Update.hh"
#include "../Aux/fileutils.hh"
#include <sstream>

namespace prf
{
    Update::Update()
    {
        isscu=true;
        isrgu=isbbu=ismcu=islocal=issequencial=false;
        st_fl=nd_fl=st_atom=nd_atom=0;
        state=0;
        itmp=0;
        nchanges=1;
        thechange.resize(1);
        nranges=natomranges=0;
        Name("A Generic Update");
        rnd=NULL;
    }

    Update::~Update() {}

    void Update::connect(Population * pl)
    {
        popl = pl;
        build_dof_list();
    }

    void Update::set(InstructionString cmd)
    {
        Logger(30)<<Name()<<"> Processing: "<<cmd.str()<<"\n";
        if (cmd.head()=="scale") setScale(strtod(cmd.tail().str().c_str(),NULL));
    }

    void Update::setScale(double gscl) {gscl;}

    void Update::configure(prf_xml::XML_Node *nd)
    {
        if (nd==NULL) return;
        for (size_t i=0;i<nd->n_children();++i) {
            prf_xml::XML_Node *dofnd=nd->child(i);
            if (dofnd==NULL or dofnd->name()!="dof") continue;
            int uid=popl->get_dof_id(dofnd->attribute("id"));
            if (uid<0) continue;
            prf_xml::XML_Node *wtnod=dofnd->child("weight");
            if (wtnod!=NULL) {
                set_weight_of_dof(uid,strtod(wtnod->value().c_str(),NULL));
            }
        }
        if (nd->child("sweep_mode")!=NULL) {
            issequencial=(nd->child("sweep_mode")->value()=="sequencial");
        }
    }

    void Update::print_setup(std::string &st)
    {
        std::ostringstream ost;
        if (issequencial) {
            ost<<(Name()+" sweeps the DOFs sequencially.\n");
        }

        for (size_t i=0;i<site_weight.size();++i) if (site_weight[i]!=1) {
            ost<<Name()<<"> using modified weight "<<site_weight[i]
                 <<" for DOF "<<dof_at_site[i].str()<<"\n";
        }
        st+=(ost.str());
    }

    size_t Update::get_site_of_dof(DOF_Info &dof)
    {
        size_t ans=0;
        for (;ans<dof_at_site.size();++ans) {
            DOF_Info sdof=dof_at_site[ans];
            if (dof.global_index==sdof.global_index or
                (dof.dof_kind==sdof.dof_kind &&
                dof.specific_global_index==sdof.specific_global_index) or
                (dof.chain==sdof.chain && dof.dof_kind==sdof.dof_kind &&
                 dof.specific_index_in_chain==sdof.specific_index_in_chain) or
                (dof.group==sdof.group && dof.dof_kind==sdof.dof_kind &&
                 dof.specific_index_in_group==sdof.specific_index_in_group))
                break;
        }
        return ans;
    }

    double Update::get_weight_of_dof(int idof)
    {
        site=get_site_of_dof(popl->get_dof_info(idof));
        if (((size_t) site)!=dof_at_site.size()) return site_weight[site];
        else return 0;
    }

    void Update::set_weight_of_dof(int idof, double vl)
    {
        size_t site=get_site_of_dof(popl->get_dof_info(idof));
        if (((size_t) site)!=dof_at_site.size()) site_weight[site]=vl;
    }

    void Update::init()
    {
        selector.set_range(0,dof_at_site.size());
        selector.set_weights(site_weight);
        selector.init();
    }

    void Update::build_dof_list() {}

    void Update::set_RandomNumberGenerator(RandomNumberBase *rn)
    {
        rnd=rn;
    }

    void Update::residue_rigid_range(int i, int &r1, int &r2)
    {
        r1=ranges[i].first;
        r2=ranges[i].second;
    }

    void Update::atom_rigid_range(int i, int &r1, int &r2)
    {
        r1=atom_ranges[i].first;
        r2=atom_ranges[i].second;
    }

    double Update::intrinsic_weight() const { return 1;}

    int Update::perform()
    {
        if (issequencial) {
            site=(site+1)%site_weight.size();
        } else site=selector.pick(rnd->shoot());
        state=1;
        //Now do something with the chosen site.
        return 1;
    }

    int Update::revert()
    {
        if (state==1) {
            for (unsigned i=0;i<nchanges;++i) {
                popl->set_dof(thechange[i].info.global_index,
                              thechange[i].before);
            }
            AtomCoordinates::revert(st_atom,nd_atom);
            state=0;
            st_fl=nd_fl=st_atom=nd_atom=0;
            return 1;
        } else {
            prf::cerr <<Name()<<"> Revert called on without the update having "
                    <<"taken place, or its having been accepted. Revert not "
                    <<"possible. \n";
            return 0;
        }
    }

    int Update::accept()
    {
        if (state==1) {
            AtomCoordinates::update(st_atom,nd_atom);
            state=0;
            st_fl=nd_fl=st_atom=nd_atom=0;
            return 1;
        } else {
            prf::cerr<<Name()<<">Accept called on without the update having "
                    <<"taken place, or its having already been accepted. "
                    <<"Accept not possible. \n";
            return 0;
        }
    }
}
