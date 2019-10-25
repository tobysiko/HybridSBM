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

#include "DihedralRestraints.hh"
#include "../../Aux/fileutils.hh"
#include "../../Aux/Constants.hh"
#include <fstream>
#include <deque>

typedef std::multimap<int,RestraintFunction *>::iterator mmiter;

DihedralRestraints::DihedralRestraints() : Energy()
{
    Name("DihedralRestraints");
    grdtyp=1;
}

DihedralRestraints::~DihedralRestraints() {}

void DihedralRestraints::set_pars(std::string s)
{
    filename=s;
}

void DihedralRestraints::init()
{
    prf_xml::XML_Node *root=prf_xml::get_xml_tree(filename);
    if (root==NULL or root->name()!="dihedral_restraints") {
        if (root==NULL) {
            prf::cerr<<Name()<<"> No valid XML tree in "<<filename<<"\n";
        } else {
            prf::cerr<<Name()<<"> Root node in "<<filename
                    <<" is not <dihedral_restraints>\n";
            delete root;
        }
        return;
    }
    root->interpret_formatted_data();
    prf::Logger blog(3);
    for (size_t i=0;i<root->n_children();++i) {
        prf_xml::XML_Node *rst=root->child(i);
        if (rst->name()!="restraint") continue;
        std::string rftype=rst->attribute("type"),dofname;
        if (rst->child("dof_id")!=NULL) dofname=rst->child("dof_id")->value();
        int uid=p->get_dof_id(dofname);
        if(uid<0) {
            prf::cerr<<Name()<<"Unknown DOF id :"<<dofname<<"\n";
            continue;
        }
        RestraintFunction *rf=NULL;
        if (rftype.empty() or rftype=="vonMises" or rftype=="v"
            or rftype=="circular_normal") {
            rf=new CircularNormal();
        } else {
            prf::cerr<<Name()<<"Unknown restraint function: "<<rftype<<"\n";
            continue;
        }
        if (rf->set_pars(rst->child("parameters"))==1) {
            c.insert(std::pair<int,RestraintFunction*>(uid,rf));
        } else {
            prf::cerr<<Name()<<"Parameter set up failed for restraint function "
                    <<rftype<<" for DOF id "<<dofname<<"\n";
            delete rf;
        }
    }
    blog<<"Created "<<c.size()<<" dihedral restraints.\n";
    if (root) delete root;
}

double DihedralRestraints::evaluate()
{
    vval=delv=0;
    for (mmiter it=c.begin();it!=c.end();++it) {
        double dof=p->get_dof(it->first);
        vval+=(*(it->second))(dof);
    }
    return vval;
}

double DihedralRestraints::deltaE(prf::Update *up)
{
    delv=0;
    for (unsigned i=0;i<up->num_changes();++i) {
        dof_change_type chng=up->change(i);
        mmiter it;
        std::pair<mmiter,mmiter> rnge;
        rnge=c.equal_range(chng.info.global_index);
        for (it=rnge.first;it!=rnge.second;++it) {
            double before,after;
            before=(*(it->second))(chng.before);
            after=(*(it->second))(chng.after);
            delv+=(after-before);
        }
    }

    return delv;
}

void DihedralRestraints::rangeEstimate(double &x1,double &x2)
{
    x1=x2=0;
    for (mmiter it=c.begin(); it!=c.end();++it) {
        x1+=(*(it->second)).estimate_min();
        x2+=(*(it->second)).estimate_max(UnivConstants::pi);
    }
}
