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

#include "ObsExpression.hh"

ObsExpression::ObsExpression() : p(NULL) {}

ObsExpression::~ObsExpression() {}

void ObsExpression::init()
{
    if (!jcobi.empty()) {
        J.set_population(p);
        J.set_dofmap(indexmap);
        gxsum.resize(3*p->NumberOfAtoms(),0.0);
        gxtmp.resize(3*p->NumberOfAtoms(),0.0);
    }
    gtmp.resize(nd,0.0);
}

void ObsExpression::reset_dof_list()
{
    indexmap.clear();
    nd=0;
}

void ObsExpression::reset_obs_list()
{
    obs.clear();
    c.clear();
    basic.clear();
    dirct.clear();
    fdelt.clear();
    jcobi.clear();
}

void ObsExpression::add_obs(double scl, Observable *o)
{
    if (o!=NULL && scl!=0) {
        size_t indx=c.size();
        if (o->grad_eval_type()==3) {
            dirct.push_back(indx);
        } else if (o->grad_eval_type()==2) {
            jcobi.push_back(indx);
        } else if (o->grad_eval_type()==1) {
            fdelt.push_back(indx);
        } else {
            basic.push_back(indx);
        }
        obs.push_back(o);
        c.push_back(scl);
    }
}

void ObsExpression::add_dof(DOF_Info &d)
{
    if (d.global_index>=0) {
        if (std::find(indexmap.begin(),indexmap.end(),d.global_index)==
            indexmap.end()) {
            indexmap.push_back(d.global_index);
            ++nd;
        }
    }
}

void ObsExpression::get_position(std::valarray<double> &x)
{
    for (size_t i=0;i<indexmap.size();++i) x[i]=p->get_dof(indexmap[i]);
}

void ObsExpression::set_position(std::valarray<double> &x)
{
    for (size_t i=0;i<x.size();++i) {
        p->set_dof(indexmap[i],x[i]);
    }
    p->Reconstruct();
//    p->EnforceBC();
    AtomCoordinates::update(0,p->NumberOfAtoms());
}

double ObsExpression::get_coordinate(int i)
{
    return p->get_dof(indexmap[i]);
}

void ObsExpression::set_coordinate(int i, double x)
{
    p->set_dof(indexmap[i],x);
    p->Reconstruct();
    p->EnforceBC();
    AtomCoordinates::update(0,p->NumberOfAtoms());
}

double ObsExpression::value()
{
    double ans=0;
    for (size_t i=0;i<obs.size();++i) {
        ans+=c[i]*obs[i]->evaluate();
    }
    return ans;
}

void ObsExpression::gradient(std::valarray<double> &g)
{
    if (!jcobi.empty()) {
        gxsum=0;
        gxtmp=0;
        J.refresh();
    }
    g=0;
    gtmp=0;

    for (size_t i=0;i<dirct.size();++i) {
        obs[dirct[i]]->gradientDOF(gtmp,indexmap);
        gtmp*=c[dirct[i]];
        g+=gtmp;
    }
    if (!jcobi.empty()) {
        for (size_t i=0;i<jcobi.size();++i) {
            obs[jcobi[i]]->gradientXYZ(gxtmp);
            gxtmp*=c[jcobi[i]];
            gxsum+=gxtmp;
        }
        J.dot(gxsum,gtmp);
        g+=gtmp;
    }
    if (!(fdelt.empty() and basic.empty())) {
        for (size_t i=0;i<fdelt.size();++i) {
            obs[fdelt[i]]->evaluate();
        }
        NDFunctional::gradient(gtmp);
        g+=gtmp;
    }
}

double ObsExpression::delta(int idir, double h)
{
    double delt=0;
    double vp=0,vm=0;
    update.pick_dof(indexmap[idir]);
    update.set_scale(h);
    update.perform();
    for (size_t i=0;i<basic.size();++i) {
        vp+=c[basic[i]]*obs[basic[i]]->evaluate();
    }
    for (size_t i=0;i<fdelt.size();++i) {
        delt+=c[fdelt[i]]*obs[fdelt[i]]->delta(&update);
    }
    update.revert();
    update.set_scale(-h);
    update.perform();
    for (size_t i=0;i<basic.size();++i) {
        vm+=c[basic[i]]*obs[basic[i]]->evaluate();
    }
    for (size_t i=0;i<fdelt.size();++i) {
        delt-=c[fdelt[i]]*obs[fdelt[i]]->delta(&update);
    }
    update.revert();
    return vp-vm+delt;
}
