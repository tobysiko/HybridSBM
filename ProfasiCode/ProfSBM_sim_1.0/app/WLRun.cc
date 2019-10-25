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
#include "WLRun.hh"
#include <Aux/Timer.hh>

using namespace prf;
using namespace UnivConstants;

WLRun::WLRun() : BasicMCRun()
{
    mc=&wl;
    progname="WLRun";
    progheader="Wang Landau iterations";

    string disable_list[]={"num_therm_cycles","use_relaxation_cycles"};

    for (int i=0;i<2;++i) optn.disable(disable_list[i]);
    maxfupdates=10;
    ncyc_per_T_updt=1;
    optn.option("n_stages","ns",1,"(number of times the WL multiplicative factor is updated)");
}

WLRun::~WLRun() {}

int WLRun::update_T()
{
    int oldTindex=curTindex;
    if (wl.check_bins() && (icyc+1)%ncyc_per_T_updt==0) {
        curTindex=wl.update_g_factor();
    }
    if (oldTindex!=curTindex) {
        H.writeRTSnapshot();
        PH.population()->SaveSnapshot(snapshot_format,
                                      pdbfil,icyc,oldTindex, esum.Value());
    }

    if (curTindex==(maxfupdates-1)) H.enableStatistics();
    if (curTindex>=maxfupdates) {
        prf::cout<<"Maximum number of updates to the f parameter "<<maxfupdates
                <<" has been reached. Wang Landau iterations have converged. \n";
        prf::cout<<"Cycle counts "<<icyc<<" to "<<MCCYC
                <<" will not be executed\n";
        icyc=MCCYC;
    }
    return curTindex;
}

int WLRun::init_MC()
{
    swtch["thermalisation"]=true;
    swtch["preliminary_relaxation"]=true; //these two must be true in WL
    H.disableStatistics();
    NTHERM=MCCYC; //to ensure that BasicMCRun does not undo the above
    wl.set_n_temps(maxfupdates);
    NTMP=maxfupdates;
    wl.init();
    return mc->Setup();
}

void WLRun::run_relaxation_cycles()
{
    prf::clog<<"Performing high temperature Monte Carlo cycles "
    <<"to get rid of obvious steric clashes. Wang Landau iterations "
    <<"can only start if the energy of the starting system is within "
    <<"the preset range\n";
    wl.bring_to_range();
    prf::cout<<"Finished relaxation cycles. Energy is now in range. "
            <<"Wang Landau iterations can begin.\n";
}

void WLRun::writeConf()
{
    BasicMCRun::writeConf();
    write_MC_averages();
}

void WLRun::write_MC_averages()
{
    char wlstats[20];
    sprintf(wlstats,"/wlstats_%d",curTindex);
    std::string wlfile=mydir+wlstats;
    wl.save_state(wlfile);
}

int WLRun::parseCommand(InstructionString s)
{
    if (s.head()=="n_stages") maxfupdates=atoi(s.tail().str().c_str());
    return BasicMCRun::parseCommand(s);
}

void WLRun::show_basic_help()
{
    prf::cout<<"WLRun: Wang Landau iterations with ProFASi.\n"
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
            <<"and the Wang Landau algorithm.\n";
    prf::cout<<"For example, you could write \n\n"
            <<"WLRun --help Population\n\n"
            <<"to see command line options about setting up the population.\n\n";
}
#ifndef PARALLEL
int main(int argc, char *argv[])
{

    WLRun interface;

    if (interface.init(argc,argv)) interface.Run();

    return 0;
}

#else
int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    WLRun interface;
    interface.set_rank_nruns(MPI::COMM_WORLD.Get_rank(),
                            MPI::COMM_WORLD.Get_size());

    if (interface.init(argc,argv)) interface.Run();

    MPI::Finalize();

    return 0;
}

#endif
#define MAINFUNC
#include "BasicMCRun.cc"
