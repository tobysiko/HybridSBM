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

#include "ObsEnergy.hh"

using namespace prf;

ObsEnergy::ObsEnergy() : Energy()
{
    Name("ObsEnergy");
    obs_name="Unknown";
    sclf=1.0;
    obs=NULL;
}

ObsEnergy::~ObsEnergy() {}

void ObsEnergy::setScaleFactor(double x) {sclf=x;}

void ObsEnergy::connectObs(Observable* o)
{
    obs=o;
    grdtyp=obs->grad_eval_type();
}

double ObsEnergy::evaluate()
{
    if (obs) {
        obs->refresh();
        vval=obs->Value();
    } else vval=0;
    vval*=sclf;
    delv=0;
    return vval;
}

double ObsEnergy::gradientXYZ(std::valarray<double> &ans)
{
    if (obs) {
        vval=obs->gradientXYZ(ans);
    } else {
        ans=0;
        vval=0;
    }
    vval*=sclf;
    ans*=sclf;
    return vval;
}

double ObsEnergy::gradientDOF(std::valarray<double> &ans,
                              std::vector<int> &indxs)
{
    if (obs) {
        vval=obs->gradientDOF(ans,indxs);
    } else {
        ans=0;
        vval=0;
    }
    vval*=sclf;
    ans*=sclf;
    return vval;
}

std::string ObsEnergy::obsName() const {return obs_name;}

void ObsEnergy::setObsName(std::string s)
{
    obs_name = s;
    Name("E_"+s);
}

void ObsEnergy::rangeEstimate(double &x1, double &x2)
{
    if (obs) obs->rangeEstimate(x1,x2);
    else {
        x1=-1;
        x2=1;
    }
    x1*=sclf;x2*=sclf;
}

void ObsEnergy::set_pars(std::string pars)
{
    std::deque<std::string> parts;
    prf_utils::split_str(pars,',',parts);
    for (size_t i=0;i<parts.size();++i) {
        std::deque<std::string> hs;
        prf_utils::split_str(parts[i],'=',hs,2);
        if (hs.size()!=2) continue;
        if (hs[0]=="scale") setScaleFactor(strtod(hs[1].c_str(),NULL));
        else if (hs[0]=="obs" or hs[0]=="watch") setObsName(hs[1]);
    }
}
