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

#include <mpi.h>
#include "ParTempRun.hh"

using namespace prf;
using namespace UnivConstants;

ParTempRun::ParTempRun() : BasicMCRun()
{
    progname="ParTempRun";
    progheader="Parallel tempering";
    gpstfil="partemp.stat";
    swtch["keep_T_history"]=false;
    mc=&ptmp;
    optn.enable("T_update_interval");
    optn.new_switch("keep_T_history","k",false,
                "(Keep a log of all accepted temperature moves. Not default)");
    Th_cache.clear();
}
/**
  \page partemprun_opts Options for temperature history
  \li \b --T_update_interval or \b -ntup : Interval in MC sweeps at which
  temperature swaps are attempted. Example : --T_update_interval 10. Note that
if it is set to 50, it means an exchange between temperatures t3 and t4
<i>every hundred sweeps</i>. This is because t3 can exchange with t4 and t2,
but not at the same time. So, a given neighbouring pair is used for every
second attempt to exchange temperatures.
  \li \b --keep_T_history or \b -k : It's a switch, and it is off by default.
  If turned on, a complete record of all accepted temperature moves is saved.
  */

ParTempRun::~ParTempRun()
{
    if (suspstate) delete [] suspstate;
}

int ParTempRun::update_T()
{
    int oldt=curTindex;
    ptmp.SwitchTemp();
    curTindex=ptmp.CurTempIndex();

    if (swtch["keep_T_history"] && curTindex!=oldt)
        Th_cache.push_back(std::make_pair(icyc,curTindex));

    return curTindex;
}

bool ParTempRun::should_suspend(bool myrec)
{
    MPI::COMM_WORLD.Allgather(&myrec,1,MPI_LOGICAL,suspstate,1,MPI_LOGICAL);

    for (int i=0;i<nruns;++i) if (suspstate[i]) return true;

    return false;
}

void ParTempRun::flush_T_data()
{
    if (!swtch["keep_T_history"]) return;

    Output tdat;

    tdat.open((mydir+"/T_history").c_str(),"a");

    for (size_t i=0;i<Th_cache.size();++i) {
        tdat<<Th_cache[i].first<<"  "<<Th_cache[i].second<<"\n";
    }

    Th_cache.clear();

    tdat.close();
}

int ParTempRun::get_span(prf_traj::Trajectory &traj,
                         unsigned long &i1, unsigned long &i2)
{
    int nerr=BasicMCRun::get_span(traj,i1,i2);
    if (nerr!=0) return 0;
    unsigned long cycmin=i1,cycmax=i2;
    MPI::COMM_WORLD.Allreduce(&i1,&cycmin,1,MPI_UNSIGNED_LONG,MPI_MAX);
    MPI::COMM_WORLD.Allreduce(&i2,&cycmax,1,MPI_UNSIGNED_LONG,MPI_MIN);
    if (cycmin!=i1 or cycmax!=i2) {
        prf::cerr<<"ParTempRun> Available cycle span for rank "<<myrank
                <<" changed from ("<<i1<<","<<i2<<") to ("
                <<cycmin<<","<<cycmax<<")\n";
        i1=cycmin;
        i2=cycmax;
    }
    return nerr;
}

void ParTempRun::set_index(int indx, unsigned long icycl)
{
    ptmp.SetTemp(indx);
    ptmp.Synchronize(icycl/ncyc_per_T_updt);
}

void ParTempRun::write_MC_averages()
{
    BasicMCRun::write_MC_averages();
    ptmp.writeTempStat(gpstfil);
    flush_T_data();
}

void ParTempRun::writeTemperatures()
{
    ptmp.writeTemperatures(mydir+"/temperature.info");
}

int ParTempRun::init_MC()
{
    int nerr=ptmp.Setup();
    suspstate=new bool[nruns];
    Th_cache.push_back(std::make_pair((unsigned long) 0,
                                      curTindex=ptmp.CurTempIndex()));
    NTMP=ptmp.number_of_temperatures();
    return nerr;
}

void ParTempRun::init_filenames()
{
    BasicMCRun::init_filenames();
    gpstfil=mydir+string("/")+gpstfil;
}

void ParTempRun::run_relaxation_cycles()
{
    size_t prevT=ptmp.CurTempIndex();
    mc->SetBeta(1.0/Trelax);   //high temperature relaxation cycles
    int ncycl=100*PH.population()->NumberOfChains();
    Logger(10)<<"Relaxing the system for "<<ncycl<<" Monte Carlo cycles "
            <<" at temperature "<<Trelax*pru_in_kelvin
            <<" Kelvins to get rid of obvious steric clashes. \n";

    for (int icyc=0;icyc<ncycl;++icyc) {
        mc->RunCycle();
        ffh.interaction_potential()->reset_total_silently();
    }

    ptmp.SetTemp(curTindex=prevT);

    Logger(10)<<"Equilibrating at temperature T_"<<prevT<<"\n";

    for (int icyc=0;icyc<ncycl;++icyc) {
        mc->RunCycle();
        ffh.interaction_potential()->reset_total_silently();
    }

    Logger(10)<<"Relaxation and equilibration "
            <<"finished before entering the main loop, and hence does "
            <<"not enter any measurements.\n";
}

int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    ParTempRun interface;
    interface.disable("rank");
    interface.disable("nruns");
    interface.set_rank_nruns(MPI::COMM_WORLD.Get_rank(),
                            MPI::COMM_WORLD.Get_size());

    if (interface.init(argc,argv)) interface.Run();

    MPI::Finalize();

    return 0;
}

#define MAINFUNC
#include "BasicMCRun.cc"

