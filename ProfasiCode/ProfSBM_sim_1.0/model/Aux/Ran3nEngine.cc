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

#include "Ran3nEngine.hh"
#include "profasi_io.hh"
#include <sstream>

using namespace prf;
const long Ran3nEngine::MBIG=1000000000;
const long Ran3nEngine::MZ=0;
const long Ran3nEngine::MSEED=161803398;
int Ran3nEngine::inext=0;
int Ran3nEngine::inextp=0;
long Ran3nEngine::ma[55]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0,0,0
                         };
double Ran3nEngine::FAC=(1.0/MBIG);
Ran3nEngine::Ran3nEngine() {Name("Ran3n");}

Ran3nEngine::~Ran3nEngine() {}

double Ran3nEngine::shoot()
{
    long mj;
    double ret_val;

    if (++inext == 55) inext=0;

    if (++inextp == 55) inextp=0;

    mj=ma[inext]-ma[inextp];

    if (mj < MZ) mj += MBIG;

    ma[inext]=mj;

    ret_val = mj*FAC;

    if (mj == 0) ret_val = FAC;

    ++ncalls;

    return ret_val;
}

void Ran3nEngine::ResetDefaultState(long seedd)
{
    long mj,mk;
    int i,ii,k;
    Logger()(5)<<Name()<<"> Setting random number seed to "<<seedd<<"\n";;
    mj=MSEED-(seedd< 0 ? -seedd : seedd);
    mj %= MBIG;
    ma[54]=mj;
    mk=1;

    for (i=1;i<=54;++i) {
        ii=(21*i) % 55;
        ma[ii-1]=mk;
        mk=mj-mk;

        if (mk < MZ) mk += MBIG;

        mj=ma[ii-1];
    }

    for (k=0;k<4;++k)
        for (i=0;i<55;++i) {
            ma[i] -= ma[(i+31) % 55];

            if (ma[i] < MZ) ma[i] += MBIG;
        }

    inext=0;

    ncalls=0;
    inextp=31;
}


