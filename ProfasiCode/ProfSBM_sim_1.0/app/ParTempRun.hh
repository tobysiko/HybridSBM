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

#ifndef ParTempRun_HH
#define ParTempRun_HH
#include <Algorithms/ParTemp.hh>
#include "BasicMCRun.hh"

using namespace prf;

/**
  \page partemp_progref Parallel tempering with PROFASI
Parallel tempering or replica exchange method is a very popular Monte Carlo
algorithm. A set of N replicas are simulated simultaneously at N different
temperatures and configurations are exchanged between different replicas with a
Metropolis like probability.

Parallel tempering is implemented in PROFASI using MPI. The executable is
called \e ParTempRun.mex, and there is no serial counterpart to this program.
Typically, you would set up a population (chains, box size etc.) using a
settings file or the command line. The procedure to do these things is exactly
the same as in canonical MC simulations (\ref bmcrun_progref), simulated
annealing (\ref siman_progref) and simulated tempering (\ref simtemp_progref).

\section options Options
Apart from all set up options (which make sense) of described in connection
with canonical MC run (\ref bmcrun_progref), parallel tempering accepts
\li \subpage partemprun_opts
\li \ref gmc_opts
\li \subpage partemp_opts


\section setup_a_partemprun What do I need for a parallel tempering run ?
\li <b>A settings file</b>: Something to tell the program what system you
want to simulate. Look up what should be in the settings file in the reference
for canonical MC runs (\ref bmcrun_progref). Also, there is an example
settings file below.
\li <b>A temperature file</b> if you are using your own temperature set. You can
use commands like "max_temperature", "num_temperatures" in the settings file
or the command line to set up a geometric series of temperatures. For any
other temperature set, you have to create a file and tell ParTempRun.mex about
it using the command "tfile".
\li <b>Files required because of your settings</b>: If you are importing a
sequence of structure from a PDB, or XML (ProFASi style) file in the settings
file, that file needs to be present.

\section example_setting Example set up
Create an empty directory and a file called "settings.cnf" in it. Put the
following lines in it.

\verbatim
log_level 10
add_chain_pdb 1 1GB1.pdb::A,41,56
box_length 100
num_cycles 5000000
cycle_length -1
tfile T.dat
conf_write_freq 1000
rt_write_freq 1000
avg_write_freq 10000
num_therm_cycles 30000
new_obs Rg rg
new_obs ProteinRMSD hvrmsd using +HV ; struc1 1GB1.pdb:1:A,41,56 ; struc2 $::A
new_obs NativenessQ Q structure 1GB1.pdb::A,41,56
\endverbatim

This sets up a simulation of the C-terminal hairpin of 1GB1.pdb. You will need
the pdb file 1GB1.pdb for the sequence and RMSD measurements. Get it from the
protein data bank (http://www.rcsb.org). Now create a temperature file called
"T.dat" with the following lines.

\verbatim
#temperature Kelvin
370
355
340
325
310
295
280
270
\endverbatim
The first line tells the program that temperatures are given in Kelvin units.
You can also write "#temperature inverted". If this line is omitted, PROFASI
internal temperature units will be assumed.

Run the program using your system's implementation of mpirun or mpiexec.
Something like:
\verbatim
mpirun -np 8 $profasi/app/bin/ParTempRun.mex
\endverbatim

If you have no idea what kind of temperatures to use, a good start is to use
a geometric series. To do this, remove the line "tfile T.dat" in the settings
file. In its place, add the following 3 lines:
\verbatim
tmin 274 Kelvin
tmax 374 Kelvin
ntmp 8
\endverbatim

If you now run the program as above, it will automatically create a simple
temperature set (in a geometric series) and use it. The same can be done by
deleting the tfile line in the settings file and instead of adding the
3 lines above, passing them in the command line:

\verbatim
mpirun -np 8 $profasi/app/bin/ParTempRun.mex -tmin 274 -tmax 374 -ntmp 8
\endverbatim

If you try to run the above example set up with 5 computing cores, it will be
an error. Each replica needs to be a separate MPI process. If you run it on 16
(or any integral multiple of 8) cores on the other hand, it will not be an
error. There will then be 2 (or more) replicas at each temperature, the
section \ref multiplexing discusses this further.

\section multiplexing Multiplexing
If you set up a parallel tempering run with 8 temperatures, and run it on
32 MPI processes, there will be 4 replicas at each temperature. At the time
of temperature exchanges, a replica at a certain temperature t0 attempts to
exchange temperature with another at a neighbouring temperature t1. In the
current example, there are 4 such replicas at t0 and 4 at t1. The exchange
attempts are done between random pairings between these two groups. Without
the random pairings, the run with 32 MPI processes is equivalent to 4
independent runs with 8 MPI processes each. Because of the random couplings,
any replica might exchange temperatures with any other in the entire run. This
idea is called "multiplexing".

The benefits of multiplexing are, at present, unclear. So, it is turned off
by default, which means that if you run 32 MPI processes with 8 temperatures,
you end up with data that is equivalent to 4 independent runs with 8 processes
each. But it can be turned on easily in the settings file with the line
"multiplexing on", or on the command line,

\verbatim
mpirun -np 8 $profasi/app/bin/ParTempRun.mex --multiplexing
\endverbatim

\section Reconstructing temperature history of a single replica
A typical ProFASi simulation is about 10 million sweeps long. Therefore it
would take a lot of disk space to store the run time history at every MC
sweep. Normally, the run time history is written only once every 1000 sweeps
or so. But often one wants to know exactly how the temperature flowed among
the different replica. That information can not be reconstructed from the
configuration files or the "rt" files without re-running the simulation. For
this reason, there is an option "--keep_T_history".

If the option "--keep_T_history" is given, a file "T_history" is created in
the output directory of every replica. That file contains (time, temperature)
pairs for every successful temperature change. So, a certain portion in that
file reads:
\verbatim
...
9543   7
9577   8
9601   7
...
\endverbatim
it means that the replica was at temperature with index 7 from sweep
number 9543 through 9576, changed to temperature 8 and remained there
from sweep 9577 through 9600 before coming back to temperature 7. This
information combined from all the replicas can be used to reconstruct a
complete history of temperature updates in the run.
<div style="color:red;"> But remember: The T_history files will be huge.
Use this option, only if you really want this information! </div>

Using the option --T_update_interval (e.g. --T_update_interval 10), one can make
parallel tempering attempt temperature exchanges every 10 Monte Carlo
sweeps instead of every sweep.
*/

//! Parallel tempering class
/**
  * The ParTempRun class adapts the BasicMCRun class to perform replica
  exchange Monte Carlo simulations. The actual replica exchange or parallel
  tempering algorithm is implemented in the ParTemp class. This class sets
  up a run using BasicMCRun functions, and adds tracking of the temperature
  movement for the replicas.

  In parallel tempering, the temperature of one replica varies in the course
  of the run, and this history can not be reconstructed based on the random
  number seed, initial conformation etc of one particular replica. It is a
  function of the history in all other replicas. Therefore, to aid data
  analysis while avoiding re-running on all the replicas, a temperature
  history file, "T_history" is created <b>on demand</b>.

  For usage information, please consult \ref partemp_progref .
  \sa \ref partemp_progref , ParTemp
  */
class ParTempRun : public BasicMCRun
{
public:
    ParTempRun();
    ~ParTempRun();
protected:
    int init_MC();
    void init_filenames();
    int update_T();
    bool should_suspend(bool myrec);
    void flush_T_data();
    void write_MC_averages();
    void writeTemperatures();
    void run_relaxation_cycles();
    void set_index(int indx, unsigned long icycl);
    int get_span(prf_traj::Trajectory &traj,
                 unsigned long &i1, unsigned long &i2);
private:
    std::string gpstfil;
    ParTemp ptmp;
    std::deque<std::pair<unsigned long,int> > Th_cache;
    bool *suspstate;
};

#endif
