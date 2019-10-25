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

#include "HamRepEx.hh"
#include <fstream>
#include <algorithm>
#include "../Aux/RandomPermutation.hh"

using std::ifstream;

using namespace prf;

using namespace prf_utils;
using std::vector;
using std::string;

HamRepEx::HamRepEx() : HGMC()
{
    nswatmpt=0;
    multiplexing=false;
    set_n_temps(4);
    explicitntmps=false;
    initTmethod=ADJACENT;
    par.new_switch("multiplexing","mul",true,"(Use multiplexing. Default)");
    par.option("initial_T_assignment","initT",1,
               "(How initial temperatures for different ranks are assigned)");
    descr="Hamiltonian Replica Exchange";
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
HamRepEx::~HamRepEx()
{
    if (Erarray) delete(Erarray);

    if (Tarray) delete(Tarray);

    if (ndarray) delete(ndarray);
}



unsigned HamRepEx::SwitchLambda()
{
	bool debugHR = false;

	double lambda_low, lambda_high;

	int l,h;

	if (my_rank==0){
		l = Tarray[n_procs-1];
		h = Tarray[my_rank+1];
	}else if (my_rank==n_procs-1){
		l = Tarray[my_rank-1];
		h = Tarray[0];
	}else{
		l = Tarray[my_rank-1];
		h = Tarray[my_rank+1];
	}
	
	lambda_low  = lambdas[ l ];
	lambda_high = lambdas[ h ];

    double emsg[4];
    emsg[0] = ffh->interaction_potential()->value();
    emsg[1] = ffh->interaction_potential()->evaluateLambda( lambda_low ); // not needed ?
    emsg[2] = ffh->interaction_potential()->evaluateLambda( lambda_high );
    
    if (debugHR) std::cout<<"HamRepEx::SwitchLambda ... rank "<< my_rank<<"  ; lower("<<l<<")=(L="<<lambda_low<<",E="<<emsg[1]<<"); this("<<Tarray[my_rank]<<")=(L="<<lambdas[Tarray[my_rank]]<<",E="<<emsg[0]<<"); higher("<<h<<")=(L="<<lambda_high<<",E="<<emsg[2]<<"); beta="<<bt<<" - (before) ScaleGoPot="<<	ffh->interaction_potential()->get_scaleSBM()<<"\n";
//gather energy from all on to node 0
    
    emsg[3]=ran->shoot();

    MPI::COMM_WORLD.Gather(&emsg,4,MPI_DOUBLE,Erarray,4,MPI::DOUBLE,0);
    
    

    if (my_rank==0) {
    	if (debugHR) {
			std::cout<<"Tarray: ";
			for (unsigned i=0; i<n_procs; i++){
				std::cout<<Tarray[i]<<", ";
			}
			std::cout<<"\n";

			std::cout<<"ndarray: ";
			for (int i=0; i<n_procs; i++){
				std::cout<<ndarray[i]<<", ";
			}
			std::cout<<"\n";
    	}
    	double Pswap, Eswapped, Ecurrent, E_i_Xj, E_j_Xi, E_i_Xi, E_j_Xj, deltaE, ep;

    	if (debugHR) prf::cout<<"Exchange cycle "<<nswatmpt<<"\n";

        int exstart = nswatmpt&1;  // if cycle even number, then 0, else 1

        if (debugHR) prf::cout<<"exstart="<<exstart<<"\n";

        for (size_t ti=exstart;ti<ntmp-1;ti+=2) {
            for (size_t i=0;i<nnodesT;++i)
            	rans[i]=Erarray[4*(nnodesT*ti+i)+1];

            if (multiplexing)
            	random_permutation<size_t>(perm,rans);

            for (size_t inod=nnodesT*ti;inod<nnodesT*(ti+1);++inod) {

            	size_t jnod=nnodesT*(ti+1)+perm[inod-nnodesT*ti];

            	if (debugHR) prf::cout<<"Exchange attempt : T-positions "<<inod<<" with "<<jnod<<"  "<<" corresponding to nodes "<<ndarray[inod]<<" and "<<ndarray[jnod]<<". ti="<<ti<< "  ";
                
                E_i_Xj = Erarray[4*ndarray[inod]+2];
                E_j_Xi = Erarray[4*ndarray[jnod]];
                Eswapped = E_i_Xj + E_j_Xi;
                
                E_i_Xi = Erarray[4*ndarray[inod]];
                E_j_Xj = Erarray[4*ndarray[jnod]+2];
                Ecurrent = E_i_Xi + E_j_Xj ;

                deltaE = Ecurrent - Eswapped;
                ep = -bt * deltaE;
                Pswap = exp( ep);

                if (debugHR) prf::cout<<E_i_Xj<<", "<<E_j_Xi<<", "<<E_i_Xi<<", "<<E_j_Xj<<", "<<Ecurrent<<", "<<Eswapped<<", "<<deltaE<<", "<<ep<<"; RN="<<Erarray[4*ndarray[jnod]+3]<<"; Pswap="<<Pswap<<"\n";
                
                if ( Erarray[4*ndarray[jnod]+3] < Pswap ) {
                    std::swap(        ndarray[inod],         ndarray[jnod]  );
                    std::swap( Tarray[ndarray[inod]], Tarray[ndarray[jnod]] );
                    if (debugHR) prf::cout<<inod<<"-"<<jnod<<" swap SUCCEEDED!\n";
                } else
                	if (debugHR) prf::cout<<inod<<"-"<<jnod<<" swap failed!\n";
            }
        }
    }

    MPI::COMM_WORLD.Scatter(Tarray,1,MPI_UNSIGNED,&itmp,1,MPI_UNSIGNED,0);

    SetLambda(itmp); // do I need this?
    ffh->interaction_potential()->set_scaleSBM(lambdas[itmp]);

    if (debugHR) std::cout<<my_rank<<" - (after) ScaleGoPot="<<	ffh->interaction_potential()->get_scaleSBM()<<" - itmp="<<itmp<<"\n";
    ++nswatmpt;
    return itmp;
}

int HamRepEx::Setup()
{
    Logger blog(20);
    nswatmpt=0;
    my_rank=MPI::COMM_WORLD.Get_rank();
    n_procs=MPI::COMM_WORLD.Get_size();
    prf::cout<<"my_rank "<<my_rank<<"\n";
    prf::cout<<"n_procs "<<n_procs<<"\n";
    if (!explicitntmps) {
        set_n_temps(n_procs);
        explicitntmps=false;
    }
    blog<<"Hamiltonian Replica Exchange setup:\n";
    blog<<"number of processes = "<<n_procs<<"\n";
    blog<<"number of temperatures = "<<ntmp<<"\n";

    int nerr=0;

    if (n_procs<ntmp || n_procs%ntmp!=0) {
        prf::cerr<<"Hamiltonian Replica Exchange: number of runs = "<<n_procs
        <<" and number of lambdas = "<<ntmp<<"\n"
        <<"Hamiltonian Replica Exchange with the ProFASi HamRepEx class can be done "
        <<"only if the number of runs is an integral multiple of the"
        <<" number of lambdas.\n";
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

        Erarray = new double[4*n_procs];
        Tarray = new unsigned[n_procs];
        ndarray = new int[n_procs];

        for (size_t i=0;i<n_procs;++i) {
            Erarray[4*i]=Erarray[4*i+1]=Erarray[4*i+2]=Erarray[4*i+3]=0;
        }
        if (initTmethod==STAGGERED) init_T_multiplex_staggered();
        else if (initTmethod==SORTED) assign_T_esort();
        else init_T_multiplex_adjacent();

        SetLambda(Tarray[my_rank]);
        ffh->interaction_potential()->set_scaleSBM(lambdas[my_rank]);

//        std::cout<<"Tarray: ";
//		for (unsigned i=0; i<n_procs; i++){
//			std::cout<<Tarray[i]<<", ";
//		}
//		std::cout<<"\n";
//
//		std::cout<<"ndarray: ";
//		for (int i=0; i<n_procs; i++){
//			std::cout<<ndarray[i]<<", ";
//		}
//		std::cout<<"\n";




        blog<<"Hamiltonian Replica Exchange set up with lambda values...\n";

        for (size_t c=0;c<ntmp;++c) {
            blog<<"lambdas "<<c<<" = "<<lambdas[c]<<"\n";
        }

    } else {
        prf::cerr<<"Hamiltonian replica exchange: "<<nerr<<" errors in Setup\n";
        prf::cerr<<"Aborting execution\n";
        prf::cerr<<"n_procs "<<n_procs<<"\n";
        MPI::COMM_WORLD.Barrier();
        exit(2);
    }
    //exit(0);
    return nerr;
}

int HamRepEx::Synchronize(int iatmpts)
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

int HamRepEx::parseCommand(InstructionString s)
{
    if (s.head()=="multiplexing") {
        if (s.tail().str()=="on") enable_multiplexing();
        else disable_multiplexing();
    } else if (s.head()=="initial_T_assignment") {
        if (s.tail().str()=="staggered") initTmethod=STAGGERED;
        else if (s.tail().str()=="sorted") initTmethod=SORTED;
        else {
            if (s.tail().str()!="ADJACENT") {
                prf::cerr<<"HamRepEx> Unknown initial temperature assignment "
                        <<"method "<<s.tail().str()<<"\n";
            }
            initTmethod=ADJACENT;
        }
    }

    return HGMC::parseCommand(s);
}

void HamRepEx::print_setup()
{
    HGMC::print_setup();
    prf::cout<<"Index   Lambda\n";

    for (size_t i=0;i<ntmp;++i) {
        prf::cout<<i<<'\t'<<lambdas[i]<<'\n';
    }

    prf::cout<<"Temperature index at rank "<<my_rank<<" is "<<itmp<<"\n";
}

int HamRepEx::init_T_multiplex_adjacent()
{
    for (unsigned i=0;i<n_procs;++i) {
        Tarray[i]=i/nnodesT;
        ndarray[i]=i;
    }
    return itmp=Tarray[my_rank];
}

int HamRepEx::init_T_multiplex_staggered()
{
    for (unsigned i=0;i<n_procs;++i) {
        Tarray[i]=i%ntmp;
        int j=nnodesT*(i%ntmp)+i/ntmp;
        ndarray[j]=i;
    }
    return itmp=Tarray[my_rank];
}

int HamRepEx::assign_T_esort()
{
	double lambda_low, lambda_high;

		int l,h;

		if (my_rank==0){
			l = Tarray[n_procs-1];
			h = Tarray[my_rank+1];
		}else if (my_rank==n_procs-1){
			l = Tarray[my_rank-1];
			h = Tarray[0];
		}else{
			l = Tarray[my_rank-1];
			h = Tarray[my_rank+1];
		}

		lambda_low  = lambdas[ l ];
		lambda_high = lambdas[ h ];
	
    double emsg[4];
    emsg[0]=ffh->interaction_potential()->value();
    emsg[1]=ffh->interaction_potential()->evaluateLambda(lambda_low);
    emsg[2]=ffh->interaction_potential()->evaluateLambda(lambda_high);
    emsg[3]=0;
    ffh->interaction_potential()->set_scaleSBM(lambdas[my_rank]);

    MPI::COMM_WORLD.Gather(&emsg,4,MPI_DOUBLE,Erarray,4,MPI::DOUBLE,0);
    
    if (my_rank==0) {

    	std::vector<std::vector<double> > emap(n_procs);

    	for (unsigned i=0;i<n_procs;++i) {
        	std::vector<double> tmp(4);
            tmp[0] = Erarray[4*i];
            tmp[1] = Erarray[4*i+1];
            tmp[2] = Erarray[4*i+2];
            tmp[3] = (double) i;
        	emap[i] = tmp;
        	//std::cout<<i<<":"<<emap[i][0]<<", "<<emap[i][1]<<", "<<emap[i][2]<<", "<<emap[i][3]<<"\n";
        }


        std::sort(emap.begin(),emap.end());

        for (size_t i=0; i<emap.size();++i) {
            unsigned tindx = (ntmp-1)-i/nnodesT;
            Tarray[(int) emap[i][3]] = tindx;

            unsigned j = nnodesT*ntmp-i-1;
            ndarray[j]=(int) emap[i][3];
            //std::cout<<i<<":"<<emap[i][0]<<", "<<emap[i][1]<<", "<<emap[i][2]<<", "<<emap[i][3]<<"\n";
        }
    }

    MPI::COMM_WORLD.Bcast(Tarray, n_procs,MPI_UNSIGNED,0);
    MPI::COMM_WORLD.Bcast(ndarray,n_procs,MPI_INT,0);

    return itmp=Tarray[my_rank];
}
