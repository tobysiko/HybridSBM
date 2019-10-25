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

#include "MC.hh"
#include <cstdlib>

using namespace prf;

using namespace prf_utils;
using std::string;

using namespace UnivConstants;

MC::MC() : HandlerBase(), cyclgt(-1), bt(2.0) {
    par.option("steps_per_cycle","lcyc",1,
               "(Number of steps in a sweep (-1 for N_dof))");
    par.option("temperature", "T",1,"(in ProFASi units or Kelvin)");
    popl=NULL;
    ffh=NULL;
    descr="Metropolis-Hastings Monte Carlo";
    debug = false;
}
/**
  \page settings_mc Commands for the Monte Carlo
  \section settings_mcpar Sweep length, temperature etc.
  \li \b --steps_per_cycle or \b -lcyc : Number of elementary Monte Carlo steps
  in a Monte Carlo cycle or sweep. The default length of a sweep in ProFASi
  v. 1.5 is the number of degrees of freedom in the system.
  \li \b --temperature or \b -T : Temperature for the MC run. You can
  specify the temperature in Kelvin: -T "300 Kelvin". If no unit is given,
  it is assumed that the temperature is in ProFASi's internal units, like in:
  -T 0.5.
  */
MC::~MC() {}

void MC::Connect(Population *p)
{
    popl=p;
    uph.autoSelect(popl);
}

Update *MC::perform_update()
{
    return uph.perform_update(0);
}

int MC::SimpleStep()
{
    Update *up=perform_update();
    double cost=ffh->interaction_potential()->deltaE(up);
    double r=-log(ran->shoot()/(up->intrinsic_weight()))/bt;

    if (cost<r) {
        ffh->interaction_potential()->accept(up);
        uph.accept_update();
        return 1;
    } else {
        ffh->interaction_potential()->reject(up);
        uph.reject_update();
        return 0;
    }
}

int MC::Step()
{
    Update *up=perform_update();
    double r=-log(ran->shoot()/(up->intrinsic_weight()))/bt;
    double cost=ffh->interaction_potential()->deltaE(up,r);

    if (cost<r) {
        ffh->interaction_potential()->accept(up);
        uph.accept_update();
        return 1;
    } else {
        ffh->interaction_potential()->reject(up);
        uph.reject_update();
        return 0;
    }
}

int MC::Setup()
{
    Logger blog;
    int nerr=0;
    if (popl==NULL) {
        prf::cerr<<"MC> Error! No population object passed to MC!\n";
        ++nerr;
    }
    if (ran==NULL) {
        prf::cerr<<"MC> Error! No random number generator is known to MC!\n";
        ++nerr;
    }
    uph.RandomNumberGenerator(ran);

    if (ffh==NULL) {
        prf::cerr<<"MC> Error! No force field handler is known to MC!\n";
        ++nerr;
    } else if (ffh->interaction_potential()==NULL) {
        prf::cerr<<"MC> Error! Force field handler returned a NULL interaction "
                <<"potential. May be the force field has not been initialized "
                <<"before MC set up was attempted ?";
        ++nerr;
    } else {
        blog(10)<<"MC> force field summary...\n";
        blog<<ffh->interaction_potential()->summary();
    }
    if (nerr!=0) return nerr;
    uph.assign_probs();
    setup_cycles();
    uph.init();
    return nerr;
}

void MC::SetBeta(double be)
{
    bt=be;
    uph.setBeta(be);
}

void MC::SetCycleLength(int ii)
{
    cyclgt=ii;
}

void MC::setup_cycles()
{
    if (cyclgt<0) {
        int ntotdof=0,nrtdof=0,nbbdof=0,nrgdof=0;

        for (int i=0;i<popl->NumberOfChains();++i) {
            nrtdof+=popl->Chain(i)->numRTdof();
            nbbdof+=popl->Chain(i)->numBBdof();
            nrgdof+=3;
        }

        int uprots(0),upbbs(0),uprgs(0);

        for (size_t i=0;i<uph.num_updates();++i) {
            Update *up=uph.used_update(i);
            if (up->sidechain_update()) ++uprots;

            if (up->backbone_update()) ++upbbs;

            if (up->rigid_chain_update()) ++uprgs;
        }

        if (uprgs==0) nrgdof=0;

        if (upbbs==0) nbbdof=0;

        if (uprots==0) nrtdof=0;

        ntotdof=nrtdof+nbbdof+nrgdof;

        cyclgt=ntotdof;
    }
}

int MC::parseCommand(InstructionString s)
{
    if (s.head()=="steps_per_cycle" || s.head()=="cycle_length")
        SetCycleLength(atoi(s.tail().str().c_str()));
    else if (s.head()=="temperature") {
    	//prf::cout << "Where does the T go?\n";
        bt=1.0/make_temperature(s.tail().str());

        std::cout<<"Beta="<<bt<<" ("<<s.tail().str()<<")\n";
        SetBeta(bt);
        //exit(1);
    }
    return 1;
}

void MC::parseCommands(std::list<InstructionString> &cmds, int argc, char *argv[])
{
    HandlerBase::parseCommands(cmds,argc,argv);
    uph.parseCommands(cmds,argc,argv);
}

void MC::show_help()
{
    HandlerBase::show_help();
    prf::cout<<"\nOther options, relating to conformational updates...\n";
    uph.show_help();
}

void MC::print_setup()
{
    prf::cout<<"\n1 Monte Carlo Cycle = "<<CycleLength()
            <<" Elementary Monte Carlo Steps or Updates\n\n";

    uph.print_setup();

    prf::cout<<"Number of Temperatures = 1\n";
    prf::cout<<"Temperature fixed at "<<temperature()<<" units or "
    <<temperature_in_kelvin()<<" Kelvin.\n";
}

std::string MC::ConfSignature()
{
    std::ostringstream ost;
    ost<<"simulation_method "<<descr<<"\n";
    ost<<"ntmp 1\n";
    ost<<"temperature "<<temperature_in_kelvin()<<" Kelvin\n";
    return ost.str();
}
