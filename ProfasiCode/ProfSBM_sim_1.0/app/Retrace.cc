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

#include "Retrace.hh"
#include <Aux/Timer.hh>

using namespace prf;

using namespace UnivConstants;
using std::min;

Retrace::Retrace() : BasicMCRun()
{
    time_to_change_T=0;
    nextTindex=0;
    Tmax=0.55;Tmin=0.41;NTMP=-8;
    Tupplan="T_history.old";
    optn.option("T_update_plan","th",1,"(Most often, a T_history file generated by a parallel tempering run)");
    mc=&gmc;
}

Retrace::~Retrace() {}

int Retrace::update_T()
{
    int ret=0;
    ranh.generator()->shoot();

    if (icyc==time_to_change_T) {
        gmc.SetTemp(curTindex=nextTindex);

        if (Tupin.good()) {
            Tupin>>time_to_change_T;
            Tupin>>nextTindex;
        }

        ret=1;
    }

    return ret;
}

int Retrace::init_MC()
{
    int nerr=gmc.Setup();
    NTMP=gmc.number_of_temperatures();
    return nerr;
}

int Retrace::init_resume()
{
    int ansbsc=BasicMCRun::init_resume();
    icyc=read_point;
    int prevTindex=0;

    if (STestFile(Tupplan.c_str())==0) {
        prf::cerr<<"Unable to open old temperature update record file "<<Tupplan<<"\n";
        exit(1);
    }

    Tupin.open(Tupplan.c_str());

    Tupin>>time_to_change_T;
    Tupin>>nextTindex;
    prf::cout<<"Found program start initial index of temperature to be "<<nextTindex<<"\n";

    if (nextTindex>=NTMP) {
        prf::cout<<"Interpreting that as "<<(nextTindex=nextTindex%NTMP)<<"\n";
    }

    do {
        Tupin>>time_to_change_T;
        prevTindex=nextTindex;
        Tupin>>nextTindex;
    } while (time_to_change_T<read_point);

    gmc.SetTemp(prevTindex%NTMP);

    prf::cout<<"Simulation will start with temperature index "<<prevTindex%NTMP<<"\n";

    prf::cout<<"Corresponding to the prevailing temperature at the beginning of sweep "
            <<read_point<<"\n";

    prf::cout<<"For comparison, while restoring the configuration to sweep "
            <<read_point<<" the temperature index was set to "<<curTindex<<"\n";

    return ansbsc;
}

void Retrace::writeTemperatures()
{
    gmc.writeTemperatures(mydir+"/temperature.info");
}

void Retrace::print_MC_setup()
{
    highlight("Retrace ...");
    prf::cout<<"\n1 Monte Carlo Cycle = "<<mc->CycleLength()
            <<" Elementary Monte Carlo Steps or Updates\n\n";
    prf::cout<<"Temperature history file = "<<Tupplan<<"\n";
}

int Retrace::parseCommand(InstructionString s)
{
    if (s.head()=="T_update_plan") Tupplan=s.tail().str();

    return BasicMCRun::parseCommand(s);
}

void Retrace::show_basic_help()
{
    prf::cout<<"Retrace: Retrace the path of an old parallel tempering run.\n"
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
            <<"Retrace --help Population\n\n"
            <<"to see command line options about setting up the population.\n\n";
}

int main(int argc, char *argv[])
{

    Retrace interface;

    if (interface.init(argc,argv)) interface.Run();

    return 0;
}

#define MAINFUNC
#include "BasicMCRun.cc"