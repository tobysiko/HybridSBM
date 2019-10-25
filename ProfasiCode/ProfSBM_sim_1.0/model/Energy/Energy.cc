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

#include "Energy.hh"

using namespace std;

namespace prf
{
    Energy::Energy() : Observable()
    {
        p=NULL;
        Name("A Generic Energy Term");
        vval=delv=0.0;
        inuse=false;
        initialized=false;
        xbin0=0.5;
    }

    Energy::~Energy() {}

    void Energy::Connect(Population *pl)
    {
        p=pl;
    }

    void Energy::refresh()
    {
        if (not inuse) evaluate();
        obsval=vval;
    }

    double Energy::evaluate() {
        delv=0;
        //Somehow calculate the energy, assign it to vval
        return vval;
    }

    double Energy::deltaE(Update *updt)
    {
        // This is the trivial way to calculate delta. Normally, there
        // is a better way to calculate it based on the properties of the
        // particular energy term and the conformational move represented
        // by the update updt.
        double eold=vval;
        evaluate();
        delv=(vval-eold);
        vval=eold;
        return delv;
    }

    double Energy::delta(Update *u)
    {
        return deltaE(u);
    }

    double Energy::deltaEwithlimit(Update *updt,double emax)
    {
        return deltaE(updt);
    }

    void Energy::Accept(Update *updt)
    {
        vval += delv;
    }

    void Energy::Revert(Update *updt) {}

    void Energy::rangeEstimate(double &x1, double &x2)
    {
        x1=-50;x2=0;
        //arbitrary, should definitely be over-written by inherittance
    }

    int Energy::init_obs()
    {
        if (Observable::init_obs()==0) return 0;
        if (!inuse) init();
        return 1;
    }

    void Energy::init() {}
    void Energy::re_init()
    {
        initialized=false;
        init();
    }

    void Energy::set_pars(std::string pars)
    {

    }
}
