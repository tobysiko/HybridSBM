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
#include "SimTempRun.hh"

using namespace prf;

using namespace UnivConstants;
using std::max;

SimTempRun::SimTempRun() : BasicMCRun()
{
    progheader="Simulated tempering";
    progname="SimTempRun";
    ogpfile="gpars.out";gpstfil="gpars.stat";
    optn.enable("T_update_interval");
    mc=&stmp;
}

SimTempRun::~SimTempRun() {}

int SimTempRun::init_MC()
{
    int nerr=stmp.Setup();
    NTMP=stmp.number_of_temperatures();
    return nerr;
}

int SimTempRun::update_T()
{
    esum.refresh();
    return curTindex=stmp.SwitchTemp();
}

void SimTempRun::write_MC_averages()
{
    BasicMCRun::write_MC_averages();
    stmp.writeNewGPars(ogpfile);
    stmp.writeTempStat(gpstfil);
}

void SimTempRun::writeTemperatures()
{
    stmp.writeTemperatures(mydir+string("/temperature.info"));
}

void SimTempRun::init_filenames()
{
    BasicMCRun::init_filenames();
    ogpfile=mydir+string("/")+ogpfile;
    gpstfil=mydir+string("/")+gpstfil;
}

void SimTempRun::show_basic_help()
{
    prf::cout<<"SimTempRun: Simulated Tempering with ProFASi.\n"
            <<"The command line options for this program are divided into\n"
            <<"different categories. Run the program with the option\n"
            <<"--help category_name(s) to get help on that category. The long\n"
            <<"forms of all command line options can be used as commands in\n"
            <<"a settings file. \n\n"
            <<"Help is available for the following categories ...\n\n"
            <<"Controls: run length, "
            <<"available time etc. \n"
            <<"Population: sequence, number of chain, box length etc. \n"
            <<"Energy: Choice of force field \n"
            <<"Random: Choice of random number generator \n"
            <<"MC: Options about conformational updates, Monte Carlo / "
            <<"simulated tempering.\n";
    prf::cout<<"For example, you could write \n\n"
            <<"SimTempRun --help Population\n\n"
            <<"to see command line options about setting up the population.\n\n";
}

#ifndef PARALLEL
int main(int argc, char *argv[])
{
    SimTempRun interface;

    if (interface.init(argc,argv)) interface.Run();

    return 0;
}

#else
int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    SimTempRun interface;
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
