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

#include "HelixSegment.hh"
#include "../Aux/Constants.hh"

using namespace UnivConstants;

HelixSegment::HelixSegment()
{
    Name("HelixSegment");
    nhturnresmax=1;
    fixed_his=true;
}

HelixSegment::~HelixSegment() {}

int HelixSegment::init_obs()
{
    if (Observable::init_obs()==0) return 0;

    nhturnresmax=0;

    for (int i = 0; i < p->NumberOfChains(); ++i) {
        int jstart=0,jend=p->Chain(i)->numAminoAcids();

        if (p->Chain(i)->NtermLigand()==NULL) ++jstart;

        if (p->Chain(i)->CtermLigand()==NULL) --jend;

        nhturnresmax+=std::max((jend-jstart-2),0);
    }


    return 1;
}

double HelixSegment::evaluate()
{
    double phi=0,psi=0;
    int hstart=-2;
    bool helixcontd=false;
    double obvl = 0;

    for (int i = 0; i < p->NumberOfChains(); ++i) {
        int jstart=0,jend=p->Chain(i)->numAminoAcids();

        if (p->Chain(i)->NtermLigand()==NULL) ++jstart;

        if (p->Chain(i)->CtermLigand()==NULL) --jend;

        for (int j = jstart; j < jend; ++j) {
            phi = p->Chain(i)->RamachandranPhi(j);
            psi = p->Chain(i)->RamachandranPsi(j);

            if ((phi<helix_phimin)||(phi>helix_phimax)
                ||(psi<helix_psimin)||(psi>helix_psimax)) {
                helixcontd=false;
            } else if (helixcontd) {
                if ((j-hstart)>1) obvl+=1;
            } else {
                helixcontd=true;
                hstart=j;
            }
        }
    }

    return obvl/nhturnresmax;
}

void HelixSegment::rangeEstimate(double &x0,double &x1)
{
    if (nhturnresmax!=0) {
        x0=-0.5/nhturnresmax;
    } else x0=0;

    x1=1-x0;
    if (!userbinsz) xbin0=1.0/nhturnresmax;
}
