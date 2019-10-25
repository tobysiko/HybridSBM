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

#ifdef PARALLEL
#include <mpi.h>
#endif
#include "SimAnnealRun.hh"

using namespace prf;
using namespace UnivConstants;

SimAnnealRun::SimAnnealRun() : BasicMCRun()
{
    progname="SimAnnealRun";
    progheader="Simulated annealing";
    mc=&siman;
}

SimAnnealRun::~SimAnnealRun() {}

int SimAnnealRun::init_MC()
{
    int nerr=siman.Setup();
    if (nerr==0) Trelax=siman.temperature(0);
    curTindex=0;

    size_t sacyclgt=0;

    for (size_t i=0;i<siman.number_of_temperatures();++i)
        sacyclgt+=siman.n_cycles_at_T(i);

    Logger blog(3);
    blog<<"It will take "<<sacyclgt<<" MC sweeps to go from the highest "
            <<"to the lowest temperature. Only an integral number of\n"
            <<"such simulated annealing passes will be performed. The "
            <<"total number of MC cycles will be adjusted if needed.\n";
    MCCYC=sacyclgt*(MCCYC/sacyclgt);
    NTMP=siman.number_of_temperatures();
    H.unsetSwitch("histogram");
    return nerr;
}

int SimAnnealRun::update_T()
{
    if (((size_t) siman.SwitchTemp())==siman.number_of_temperatures()) {
        siman.SetTemp(curTindex=0);
        re_init();
    }
    return curTindex=siman.CurTempIndex();
}

int SimAnnealRun::re_init()
{
    PH.init_coords();
    PH.reconstruct();
    ffh.interaction_potential()->reset_total_silently();

    if (swtch[string("preliminary_relaxation")]) run_relaxation_cycles();

    return 1;
}

int SimAnnealRun::parseCommand(InstructionString s)
{
    siman.parseCommand(s);
    return BasicMCRun::parseCommand(s);
}

void SimAnnealRun::writeTemperatures()
{
    siman.writeTemperatures(mydir+string("/temperature.info"));
}

void SimAnnealRun::show_basic_help()
{
    prf::cout<<"SimAnnealRun: Simulated Annealing with ProFASi.\n"
            <<"The command line options for this program are divided into\n"
            <<"different categories. Run the program with the option\n"
            <<"--help category_name(s) to get help on that category. The long\n"
            <<"forms of all command line options can be used as commands in\n"
            <<"a settings file. \n\n"
            <<"Help is available for the following categories ...\n\n"
            <<"Controls: simulated annealing settings, run length, "
            <<"available time etc. \n"
            <<"Population: sequence, number of chain, box length etc. \n"
            <<"Energy: Choice of force field \n"
            <<"Random: Choice of random number generator \n"
            <<"MC: Options about Monte Carlo / conformational updates.\n";
    prf::cout<<"For example, you could write \n\n"
            <<"SimAnnealRun --help Population\n\n"
            <<"to see command line options about setting up the population. \n\n";
}

#ifndef PARALLEL
int main(int argc, char *argv[])
{
    SimAnnealRun interface;

    if (interface.init(argc,argv)) interface.Run();

    return 0;
}

#else
int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    SimAnnealRun interface;
    interface.disable("rank");
    interface.disable("nruns");
    interface.set_rank_nruns(MPI::COMM_WORLD.Get_rank(),
                             MPI::COMM_WORLD.Get_size());
    if (interface.init(argc,argv)) interface.Run();

    MPI::Finalize();

    return 0;
}

#endif
#define MAINFUNC
#include "BasicMCRun.cc"
