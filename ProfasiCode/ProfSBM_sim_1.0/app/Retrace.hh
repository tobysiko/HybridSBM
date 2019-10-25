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

#ifndef Retrace_HH
#define Retrace_HH
#include "BasicMCRun.hh"
#include <Algorithms/GMC.hh>

/**
\page retrace Retrace the trajectory of a single replica of a parallel tempering run

Imagine you have a parallel tempering run with 64 replicas. Something really
interesting seems to happen between MC sweeps 150 million and 165 million on
replica 19. You want to follow it more closely, by generating a more densely
recorded run-time history file, or perhaps by making a movie of that part of
the trajectory. You want snapshots spaced at shorter intervals compared to
how often the original run saved its configuration. How would you do it ?

One way is to restart the run at sweep 150 million with a new settings file
that writes out rt and conf files more frequently. But since the journey
across temperature space of a replica in a parallel tempering simulation
depends on temperatures and energies in all replicas, you would need to
re-start a 64 replica run. Can one avoid the wastage of resources by somehow
generating only the desired information, with only a re-run of the replica 19 ?

This is the purpose of the program \tt Retrace. For this to work, the parallel tempering
run must have been started with the option "keep_T_history on" option in the
settings file. With this option, a sufficiently detailed history of the
journey of each replica through different temperatures is maintained in
files n0/T_history, n1/T_history etc. For each replica, this file records
what temperature indices were assigned to it after a given number of sweeps.
There is an entry in that file for every successful temperature update. It
could look something like this...

\verbatim
0  1
2  0
4  1
6  0
8  1
9  2
10  3
11  4
12  5
...
...
\endverbatim

The first column contains MC sweeps while the second contains temperature
indices. The first row indicates the initial value of the temperature
index for the current replica. From the second row on, the second column
is the temperature index assigned to the replica at the \em end of the
sweeps in column 1. In this example, the replica was initialized to
index 0. After sweep 0, it switched to temperature index 1. Then after
sweep 2 it went back to index 0. This automatically means, no new index
was assigned after sweep 1.

There is sufficient information in this history file to reconstruct
the entire temperature history of this replica. Therefore, given this
file, one can run something like a canonical Monte Carlo with changing
temperature assignments guided by these values, and retrace the trajectory
of the single replica without running the other 63 simultaneously. That
is how this program works.

In addition to the commands described in connection with parallel tempering
(\ref partemp_progref), this program interprets one more instruction:

\li \b --T_update_plan or \b -th : This sets the parallel tempering
temperature history file to use for retracing. See the documentation of
parallel tempering (\ref partemp_progref look up the bits about the option
"keep_T_history") for details on how to tell parallel tempering to remember
its temperature history.
*/

//! A program to retrace the path of an individual replica of parallel tempering
/**
  The purpose of this class is to follow the path of a single replica of
  a parallel tempering run, by using its temperature history file for
  temperature updates. It is implemented using the GMC class as the Monte
  Carlo algorithm.
  */
class Retrace : public BasicMCRun
{
public:
    Retrace();
    virtual ~Retrace();
    //! A few commands in addition to those taken by "BasicMCRun"
    int parseCommand(InstructionString s);
    void show_basic_help();
    int init_resume();
protected:
    int init_MC();
    int update_T();
    void writeTemperatures();
    void print_MC_setup();
    std::string Tupplan, tdistfile;
    double Tmin,Tmax;
    unsigned long time_to_change_T;
    int nextTindex, ncyc_per_T_updt;
    std::vector<double> Tarray;
    std::ifstream Tupin;
    GMC gmc;
};

#endif
