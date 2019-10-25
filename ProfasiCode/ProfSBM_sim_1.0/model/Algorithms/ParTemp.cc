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

#include "ParTemp.hh"
#include <fstream>
#include <algorithm>
#include "../Aux/RandomPermutation.hh"

using std::ifstream;

using namespace prf;

using namespace prf_utils;
using std::vector;
using std::string;

ParTemp::ParTemp() : GMC()
{
    nswatmpt=0;
    multiplexing=false;
    set_n_temps(4);
    explicitntmps=false;
    initTmethod=ADJACENT;
    par.new_switch("multiplexing","mul",true,"(Use multiplexing. Default)");
    par.option("initial_T_assignment","initT",1,
               "(How initial temperatures for different ranks are assigned)");
    descr="Parallel tempering";
}
/**
  \page partemp_opts Options peculiar to parallel tempering
  \li \b --multiplexing or \b -mul : Switch on multiplexing
  \li \b --initial_T_assignment or \b -initT  : Initial temperature assignment
      in parallel tempering. The possible values are "adjacent","staggered"
      and "sorted". For 4 replicas with 2 fold multiplexing: "adjacent" means
      ranks 0, 1, 2 and 3 get temperature indices 0,0,1,1. "staggered" means
      ranks 0, 1, 2 and 3 get temperature indices 0,1,0,1. "sorted" means the
      temperature indices are decided by sorting the energies of the ranks, and
      giving the lowest temperatures (highest indices) to the lowest energies.
  */
ParTemp::~ParTemp()
{
    if (Erarray) delete(Erarray);

    if (Tarray) delete(Tarray);

    if (ndarray) delete(ndarray);
}

unsigned ParTemp::SwitchTemp()
{
//gather energy from all on to node 0
    double emsg[2];
    emsg[0]=ffh->interaction_potential()->value();
    emsg[1]=ran->shoot();
    MPI::COMM_WORLD.Gather(&emsg,2,MPI_DOUBLE,Erarray,2,MPI::DOUBLE,0);

    if (my_rank==0) {
//        prf::cout<<"Exchange cycle "<<nswatmpt<<"\n";
        int exstart=nswatmpt&1;

        for (size_t ti=exstart;ti<ntmp-1;ti+=2) {
            for (size_t i=0;i<nnodesT;++i) rans[i]=Erarray[2*(nnodesT*ti+i)+1];

            if (multiplexing) random_permutation<size_t>(perm,rans);

            for (size_t inod=nnodesT*ti;inod<nnodesT*(ti+1);++inod) {
                size_t jnod=nnodesT*(ti+1)+perm[inod-nnodesT*ti];
//                prf::cout<<"Exchange attempt : T-positions "<<inod<<" with "<<jnod<<"  "
//                <<" corresponding to nodes "<<ndarray[inod]<<" and "<<ndarray[jnod]<<"  ";

                if (Erarray[2*ndarray[jnod]+1] <
                    exp((Erarray[2*ndarray[inod]]-Erarray[2*ndarray[jnod]])*
                        (beta[ti]-beta[ti+1]))) {
                    std::swap(ndarray[inod],ndarray[jnod]);
                    std::swap(Tarray[ndarray[inod]],Tarray[ndarray[jnod]]);
//                    prf::cout<<"succeeded!\n";
                } //else prf::cout<<"failed!\n";
            }
        }
    }

    MPI::COMM_WORLD.Scatter(Tarray,1,MPI_UNSIGNED,&itmp,1,MPI_UNSIGNED,0);

    SetTemp(itmp);
    ++nswatmpt;
    return itmp;
}



int ParTemp::Setup()
{
    Logger blog(20);
    nswatmpt=0;
    my_rank=MPI::COMM_WORLD.Get_rank();
    n_procs=MPI::COMM_WORLD.Get_size();
    prf::cout<<"n_procs "<<n_procs<<"\n";
    if (!explicitntmps) {
        set_n_temps(n_procs);
        explicitntmps=false;
    }
    blog<<"Parallel tempering setup:\n";
    blog<<"number of processes = "<<n_procs<<"\n";
    blog<<"number of temperatures = "<<ntmp<<"\n";

    int nerr=0;

    if (n_procs<ntmp || n_procs%ntmp!=0) {
        prf::cerr<<"Parallel Tempering: number of runs = "<<n_procs
        <<" and number of temperatures = "<<ntmp<<"\n"
        <<"Parallel tempering with the ProFASi ParTemp class can be done "
        <<"only if the number of runs is an integral multiple of the"
        <<" number of temperatures.\n";
        ++nerr;
    }

    MC::Setup();
    // MC set up in parallel tempering must be delayed until number of
    // temperature is decided.

    uph.assign_probs();

    if (nerr==0) {
        nnodesT=n_procs/ntmp;

        if (nnodesT==1) multiplexing=false;

        perm.resize(nnodesT,0);
        rans.resize(nnodesT,0);

        for (size_t i=0;i<nnodesT;++i) perm[i]=i;

        Erarray = new double[2*n_procs];
        Tarray = new unsigned[n_procs];
        ndarray = new int[n_procs];

        for (size_t i=0;i<n_procs;++i) {
            Erarray[2*i]=Erarray[2*i+1]=0;
        }
        if (initTmethod==STAGGERED) init_T_multiplex_staggered();
        else if (initTmethod==SORTED) assign_T_esort();
        else init_T_multiplex_adjacent();

        SetTemp(Tarray[my_rank]);
        blog<<"Parallel tempering set up with temperatures...\n";

        for (size_t c=0;c<ntmp;++c) {
            blog<<"beta "<<c<<" = "<<beta[c]<<"\n";
        }

    } else {
        prf::cerr<<"Parallel Tempering: "<<nerr<<" errors in Setup\n";
        prf::cerr<<"Aborting execution\n";
        prf::cerr<<"n_procs "<<n_procs<<"\n";
        MPI::COMM_WORLD.Barrier();
        exit(2);
    }
    return nerr;
}

int ParTemp::Synchronize(int iatmpts)
{
    nswatmpt=iatmpts;
    MPI::COMM_WORLD.Gather(&itmp,1,MPI_INT,Tarray,1,MPI_INT,0);
    size_t igr=0;

    for (size_t i=0;i<ntmp;++i) {
        for (size_t j=0;j<n_procs;++j) {
            if (Tarray[j]==i) ndarray[igr++]=j;

            if (igr==((i+1)*nnodesT)) break;
        }
    }

    return 1;
}

int ParTemp::parseCommand(InstructionString s)
{
    if (s.head()=="multiplexing") {
        if (s.tail().str()=="on") enable_multiplexing();
        else disable_multiplexing();
    } else if (s.head()=="initial_T_assignment") {
        if (s.tail().str()=="staggered") initTmethod=STAGGERED;
        else if (s.tail().str()=="sorted") initTmethod=SORTED;
        else {
            if (s.tail().str()!="ADJACENT") {
                prf::cerr<<"ParTemp> Unknown initial temperature assignment "
                        <<"method "<<s.tail().str()<<"\n";
            }
            initTmethod=ADJACENT;
        }
    }

    return GMC::parseCommand(s);
}

void ParTemp::print_setup()
{
    GMC::print_setup();
    prf::cout<<"Index   T.(model)\tT.(kelvin)\n";

    for (size_t i=0;i<ntmp;++i) {
        prf::cout<<i<<'\t'<<1.0/beta[i]<<'\t'
        <<UnivConstants::pru_in_kelvin/beta[i]<<'\n';
    }

    prf::cout<<"Temperature index at rank "<<my_rank<<" is "<<itmp<<"\n";
}

int ParTemp::init_T_multiplex_adjacent()
{
    for (unsigned i=0;i<n_procs;++i) {
        Tarray[i]=i/nnodesT;
        ndarray[i]=i;
    }
    return itmp=Tarray[my_rank];
}

int ParTemp::init_T_multiplex_staggered()
{
    for (unsigned i=0;i<n_procs;++i) {
        Tarray[i]=i%ntmp;
        int j=nnodesT*(i%ntmp)+i/ntmp;
        ndarray[j]=i;
    }
    return itmp=Tarray[my_rank];
}

int ParTemp::assign_T_esort()
{
    double emsg[2];
    emsg[0]=ffh->interaction_potential()->value();
    emsg[1]=0;
    MPI::COMM_WORLD.Gather(&emsg,2,MPI_DOUBLE,Erarray,2,MPI::DOUBLE,0);
    if (my_rank==0) {
        std::vector<std::pair<double, unsigned> > emap(n_procs);
        for (unsigned i=0;i<n_procs;++i) {
            emap[i].first=Erarray[2*i];
            emap[i].second=i;
        }
        std::sort(emap.begin(),emap.end());
        for (size_t i=0;i<emap.size();++i) {
            unsigned tindx=(ntmp-1)-i/nnodesT;
            Tarray[emap[i].second]=tindx;
            unsigned j=nnodesT*ntmp-i-1;
            ndarray[j]=(int) emap[i].second;
        }
    }
    MPI::COMM_WORLD.Bcast(Tarray,n_procs,MPI_UNSIGNED,0);
    MPI::COMM_WORLD.Bcast(ndarray,n_procs,MPI_INT,0);
    return itmp=Tarray[my_rank];
}
