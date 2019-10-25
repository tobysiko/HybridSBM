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

#include "SimTemp.hh"

using namespace prf;

SimTemp::SimTemp() : nswatmpt(0)
{
    par.option("g_parameter_file","gf",1,"(simulated tempering parameters)");
    gpfile="gpars.in";
    descr="Simulated tempering";
}

SimTemp::~SimTemp() {}

unsigned SimTemp::SwitchTemp()
{
    double e=ffh->interaction_potential()->value();
    unsigned nind;
    ++nswatmpt;
    nind=(ran->shootBit()) ? itmp+1:itmp-1;

    if (nind<ntmp) {
        if (ran->shoot()<exp(-((beta[nind]-beta[itmp])*e+g[nind]-g[itmp]))) {
            itmp=nind;
        }
    }
    SetTemp(itmp);

    return itmp;
}

void SimTemp::writeNewGPars(std::string flnm)
{
    Output ost(flnm.c_str());
    double gl=g[ntmp-1];
    gl+=log(std::max(((double)(dndT[ntmp-1]))/nswatmpt,0.01/ntmp));

    for (size_t i=0;i<ntmp;++i) {
        double gtmp=g[i];
        double gnew=gtmp;

        if (nswatmpt!=0) gnew+=log(std::max(((double)(dndT[i]))/nswatmpt,
                                           0.01/ntmp))-gl;

        ost<<i<<"   "<<gnew<<"\n";
    }

    ost.close();
}

void SimTemp::read_g_pars(std::string gfl)
{
    Logger blog(10);
    g.resize(ntmp,0);
    if (prf_utils::TestFile_r(gfl.c_str())==0) {
        prf::cerr<<"Simulated tempering> No g-parameter file! All g values "
        <<"set to 0. This run will probably not visit many temperatures.\n";
    } else {
        int idummy;
        std::ifstream fin(gpfile.c_str());

        for (size_t i=0;i<ntmp;++i) {fin>>idummy;fin>>g[i];}

        fin.close();

        blog <<"Read in g values from file "<<gpfile<<"\n";
    }

    blog <<"index : inv. temperature : g \n";
    for (size_t i=0;i<ntmp;++i)
        blog <<i<<"\t : "<<beta[i]<<"\t\t"<<g[i]<<"\n";
}

int SimTemp::Setup()
{
    int nerr=0;
    nerr+=MC::Setup();
    nswatmpt=itmp=0;

    SetTemp(itmp);
    if (nerr==0) read_g_pars(gpfile);

    return nerr;
}

int SimTemp::parseCommand(InstructionString s)
{
    if (s.head()=="g_parameter_file") {
        gpfile=s.tail().str();
    }
    return GMC::parseCommand(s);
}

void SimTemp::print_setup()
{
    GMC::print_setup();
    prf::cout<<"Index   T.(model)\tT.(kelvin)\tg(T) \n";

    for (size_t i=0;i<ntmp;++i) {
        prf::cout<<i<<'\t'<<1.0/beta[i]<<'\t'
        <<UnivConstants::pru_in_kelvin/beta[i]<<'\t'<<g[i]<<'\n';
    }
}

