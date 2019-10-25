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

#include "DistanceRestraints.hh"
#include <fstream>
#include <deque>
#include "../../Aux/fileutils.hh"

DistanceRestraint::DistanceRestraint() : atom1(-1),atom2(-1),f(NULL) {}

// f is deliberately not destroyed in the destructor here!
DistanceRestraint::~DistanceRestraint() {}

DistanceRestraint::DistanceRestraint(const DistanceRestraint &dr)
{
    atom1=dr.atom1;
    atom2=dr.atom2;
    f=dr.f;
}

DistanceRestraint & DistanceRestraint::operator =(const DistanceRestraint &dr)
                                                 {
    if (this!=&dr) {
        atom1=dr.atom1;
        atom2=dr.atom2;
        f=dr.f;
    }
    return *this;
}

void DistanceRestraint::delete_restraint_function()
{
    if (f) delete f;
    f=NULL;
}

double DistanceRestraint::evaluate()
{
    double dist=AtomCoordinates::dist(atom1,atom2);
    return (*f)(dist);
}

double DistanceRestraint::estimate_upper_bound()
{
    double dist=0.5*AtomCoordinates::boxL();
    return (*f)(dist);
}

int DistanceRestraint::set_pars(prf_xml::XML_Node *pars,Population *p)
{
    prf::Logger blog(10);
    atom1=atom2=-1;
    if (f!=NULL) delete f;
    f=NULL;
    if (pars==NULL or pars->name()!="restraint") return 0;
    if (pars->child("atom1")!=NULL) atom1=get_aid(pars->child("atom1")->value(),p);
    else if (pars->child("atom_1")!=NULL)
        atom1=get_aid(pars->child("atom_1")->value(),p);
    if (pars->child("atom2")!=NULL) atom2=get_aid(pars->child("atom2")->value(),p);
    else if (pars->child("atom_2")!=NULL)
        atom2=get_aid(pars->child("atom_2")->value(),p);

    std::string rftype=pars->attribute("type");
    if (rftype.empty() or rftype=="0" or rftype=="quadratic") {
        f=new RestraintFunction();
    } else if (rftype=="1" or rftype=="power_law" or rftype=="powerlaw") {
        f=new PowerLawRestraint();
    } else if (rftype=="2" or rftype=="flattened_power_law") {
        f=new FlattenedPL();
    } else if (rftype=="3" or rftype=="gaussian" or rftype=="Gaussian") {
        f=new GaussianRestraint();
    }
    if (atom1>0 and atom2>0 and f!=NULL) {
        if (f->set_pars(pars->child("parameters"))==0) {
            prf::cerr<<"Error while trying to set parameters for "<<rftype
                    <<" between atoms "<<atom1<<" and "<<atom2<<"\n";
        } else {
            blog<<"Created "<<rftype<<" restraint between atoms with unique ids "
                    <<atom1<<" and "<<atom2<<"\n";
            return 1;
        }
    }
    if (f!=NULL) delete f;
    return 0;
}

int DistanceRestraint::get_aid(std::string atm,Population *p)
{
    std::deque<std::string> parts;
    prf_utils::split_str(atm,'/',parts,4);
    int chid,resid;
    chid=atoi(parts[0].c_str());
    resid=atoi(parts[1].c_str());
    if (chid<0 or chid>=p->NumberOfChains()) {
        prf::cerr<<"Invalid chain id "<<chid<<" in restraint specification"
                <<" :"<<atm<<"\n";
        return -1;
    }
    if (resid<0 or resid>=p->Chain(chid)->numLigands()) {
        prf::cerr<<"Invalid residue number "<<resid<<" in restraint "
                <<"specification :"<<atm<<"\n"
                <<"Chain "<<chid<<" has "<<p->Chain(chid)->numLigands()
                <<" groups.\n";
        return -1;
    }
    std::string aa=p->Chain(chid)->memberLigand(resid)->TLC();
    if (aa!=parts[2]) {
        prf::cerr<<"Group "<<resid<<" in chain "<<chid<<" is a "<<aa
                <<" where as it is given as a "<<parts[2]<<" in restraint "
                <<"specification :"<<atm<<"\n";
        return -1;
    }
    std::string atmname=parts[3];
    for (size_t i=0;i<atmname.size();++i) if (atmname[i]=='_') atmname[i]=' ';
    Atom *tmpatm=p->Chain(chid)->memberLigand(resid)->labeled_atom(atmname);
    if (tmpatm!=NULL) {
        return tmpatm->UniqueId();
    }
    return -1;
}

DistanceRestraints::DistanceRestraints() : Energy()
{
    Name("DistanceRestraints");
}

DistanceRestraints::~DistanceRestraints()
{
    for (size_t i=0;i<c.size();++i) {
        c[i].delete_restraint_function();
    }
}

void DistanceRestraints::set_pars(std::string s)
{
    filename=s;
}

void DistanceRestraints::init()
{
    prf_xml::XML_Node *root=prf_xml::get_xml_tree(filename);
    if (root==NULL or root->name()!="distance_restraints") {
        if (root==NULL) {
            prf::cerr<<Name()<<"> No valid XML tree in "<<filename<<"\n";
        } else {
            prf::cerr<<Name()<<"> Root node in "<<filename
                    <<" is not <distance_restraints>\n";
            delete root;
        }
        return;
    }
    root->interpret_formatted_data();
    std::deque<DistanceRestraint> tmp;
    prf::Logger blog(3);
    for (size_t i=0;i<root->n_children();++i) {
        prf_xml::XML_Node *rst=root->child(i);
        if (rst->name()!="restraint") continue;
        DistanceRestraint res;
        if (res.set_pars(rst,p)==1) tmp.push_back(res);
    }
    c.resize(tmp.size());
    c.assign(tmp.begin(),tmp.end());
    blog<<Name()<<"> Created "<<c.size()<<" distance restraints.\n";
    if (root) delete root;
}

double DistanceRestraints::evaluate()
{
    vval=delv=0;
    for (size_t i=0;i<c.size();++i) {
        vval+=c[i].evaluate();
    }
    return vval;
}

void DistanceRestraints::rangeEstimate(double &x1,double &x2)
{
    x1=x2=0;
    for (size_t i=0;i<c.size();++i) {
        x2+=c[i].estimate_upper_bound();
    }
}

