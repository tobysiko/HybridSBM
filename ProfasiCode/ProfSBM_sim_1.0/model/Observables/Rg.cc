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

#include "Rg.hh"

using namespace prf_utils;

using namespace prf;
using std::string;
using std::vector;

Rg::Rg()
{
    Name("Rg");
    p=NULL;
    ich=-1;
    grdtyp=2;
}

Rg::~Rg() {}

/**
*\page opt_Rg Radius of gyration
\section options Available options
\li \b of_chain of_chain 2<br> Limits the calculation to the third chain in the population. If radius of gyration measurements are not limited to a single chain in this manner, it is calculated over the entire system.

\sa prf::Rg
*/
int Rg::init_obs()
{
    if (Observable::init_obs()==0) return 0;

    vector<string> parts;

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);

        if (parts[0]==string("of_chain") && parts.size()>=2) {
            int tmpint=atoi(parts[1].c_str());

            if ((tmpint<0) || (tmpint > p->NumberOfChains())) tmpint=-1;
            else Logger(log_thres)<<Name()
                <<"> Restricting calculation to chain "<<tmpint<<"\n";

            of_chain(tmpint);
        }
    }

    return 1;
}

double Rg::evaluate()
{
    int nr1=0,nr2=0;

    if (ich==-1) {
        nr1=0;
        nr2=p->NumberOfAtoms();
    } else {
        nr1=p->Chain(ich)->begin_atom();
        nr2=p->Chain(ich)->end_atom();
    }

    int nterms=0;

    double r2sum=0;
    Vector3 rsum,vtmp;

    for (int i=nr1;i<nr2;++i) {
        if (p->SpeciesOf(i) !=hydrogen) {
            rsum+= (vtmp=AtomCoordinates::vec(i));
            r2sum+=vtmp.mag2();
            ++nterms;
        }
    }

    return sqrt(r2sum/nterms- ((1.0/nterms) *rsum).mag2());
}

double Rg::gradientXYZ(std::valarray<double> &ans)
{
    int nr1=0,nr2=0;

    if (ich==-1) {
        nr1=0;
        nr2=p->NumberOfAtoms();
    } else {
        nr1=p->Chain(ich)->begin_atom();
        nr2=p->Chain(ich)->end_atom();
    }

    ans=0;
    int nterms=0;

    double r2sum=0;
    Vector3 rsum,vtmp;

    for (int i=nr1;i<nr2;++i) {
        if (p->SpeciesOf(i) !=hydrogen) {
            rsum+= (vtmp=AtomCoordinates::vec(i));
            r2sum+=vtmp.mag2();
            ans[3*i]+=vtmp.x();
            ans[3*i+1]+=vtmp.y();
            ans[3*i+2]+=vtmp.z();
            ++nterms;
        }
    }
    rsum=(1.0/nterms)*rsum;
    r2sum/=nterms;
    for (int i=nr1;i<nr2;++i) {
        if (p->SpeciesOf(i) !=hydrogen) {
            ans[3*i]-=rsum.x();
            ans[3*i+1]-=rsum.y();
            ans[3*i+2]-=rsum.z();
        }
    }
    obsval= sqrt(r2sum - rsum.mag2());
    ans*=(1.0/nterms/obsval);

    return obsval;
}

void Rg::rangeEstimate(double &x1,double &x2)
{
    if (p==NULL) {
        x1=0;
        x2=100;
    } else {
        x1=0;
        x2=3.6*sqrt((double) p->NumberOfLigands());
    }
    if (!userbinsz) xbin0=0.25;
}
