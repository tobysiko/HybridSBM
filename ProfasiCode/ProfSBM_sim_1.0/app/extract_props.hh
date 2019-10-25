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

#include <Elements/PopulationHandler.hh>
#include <Energy/FFHandler.hh>
#include <Observables/ObsHandler.hh>
#include <Aux/Trajectory.hh>

using namespace prf;

using namespace UnivConstants;

//! Calculating run-time properties from traj files
/**
This program works with a traj file generated during a run to calculate
values of various observables at the points where configurations were
saved during the run. In other words, it can regenerate the rt file
if the frequency of writing into that file is the same as the frequency
of writing configuration snapshots. It is useful, if, for instance, we want
to find the value of a new measurable quantity at all the points along
a run-time history. The program simply reads the traj file, sequencially
restores the population to each of the configurations saved during the run
(within a range, if specified), and calculates any observables as specified.

PROFASI "traj" files are meta-data files containing information about how to
find configuration snapshots at different MC cycles during a run. They don't
contain the actual data on the state of the population during the run. The
data is stored in binary format in a bunch of "conf..." files in the same
directory as the traj file.

In its present form, this program has two intended uses.

\li To reconstruct a run-time history file, rt, from a traj file and instructions
in the settings file. Imagine the situation that you did not have the NativenessQ
observable in the settings file when you started a run. While analysing data,
you realize it might be interesting and want to see how it correlates with
different energy terms. You can do that if you had the variable printed in the
rt file, but since it was not in the settings file when the run started, it will
be absent from the rt file. With this program, you can reconstruct a new rt file
(give it a different name), with the additional measurement (See "Usage" below).

\li Retrieve a series of PDB snapshots in a certain range of Monte Carlo cycles.
This is mostly useful for creating animations. We warn that Monte Carlo "time"
evolution should not be carelessly interpreted as "time" evolution. But still, an
animation is a good way to visualize what was happening in a simulation.

\section syntax Syntax
extract_props [OPTIONS] traj_file \n \n
where, [OPTIONS] could be one or more of:

\li \b --get_pdbs or \b -pdb Retrieve snapshot PDB files in a specified range
\li \b --single_file or \b -sf Store snapshots from different stored
configurations as different models in a single PDB file, instead of one file per
structure.
\li \b --get_rt or \b -rt Reconstruct run-time system properties in a range
\li \b --get_averages or \b -avg Generate averages file and profile data
\li \b --output_prefix or \b -op Needs an argument specifying what should be
prefixed to the output file names. For instance, if you are generating pdb files,
they will be named frame_1.pdb, frame_2.pdb etc. But with this option, you can
give it a prefix to attach to the filenames. If you pass "-op traj3_" to the
program, the PDB filenames will be traj3_frame_1.pdb, traj3_frame_2.pdb etc.
"-op d0/" will try to write to output pdb files called d0/frame_1.pdb,
d0/frame_2.pdb etc, which means the directory d0 should exist.
\li \b --start or \b -a The cycle number to start at.
\li \b --end or \b -z The cycle number (inclusive) to end at.
\li \b --every or \b -ev Use every n'th snapshot
\li \b --number_of_temperatures or \b -nt In some situations this program
needs to know the total number of temperatures in the simulation which created
the trajectory and it can not figure it out on its own.

If you need a single PDB file at a certain point in a run, this program could be
a sledge-hammer to kill a mosquito. Use the program extract_snapshot instead.

\note Although this program reads a settings file, it only cares about
instructions about Observables in the settings file. To interpret the contents
of the binary "conf..." files (dereferenced in the "traj" file), it uses the
header section of that file or information written in a "conf.info" file in the
same directory as the conf file. In other words, whatever you write in the
settings file about number of chains and sequences has no effect on this program.

\section usage Usage
\subsection usage1 Re-generating run-time information with new Observables
Situation: you have a finished run of some system. While analysing the data you
realize, you want to analyze the correlation of some new observables which
you did not think of when you did the run. For instance, in a two chain run,
you notice all minimum energy PDB files showing the chains either both folded
or both unfolded. You want a scatter plot with the radius of gyration of each
individual chain to see how they correlate, and perhaps to filter the scatter
plot based on values of other observables. But radius of gyration of the
individual chains was not in your settings file when you started the run. You
can use extract_props to retrieve the missing information. The two lines
you would need in the settings file for the two missing Rg observables are:

\verbatim
new_obs Rg rg0 of_chain 0
new_obs Rg rg1 of_chain 1
\endverbatim

After adding these lines (or any other observable related lines), run one of the
following:
\verbatim
$ extract_props -op tmp_ -rt n0/traj
$ extract_props -op tmp_ -rt -avg n0/traj
\endverbatim
The first example above generates the run time history file ("rt") from the
trajectory file n0/traj. The newly generated rt file will have an
\em output_prefix "tmp_", i.e., it will be called "tmp_rt". A file called
"tmp_rtkey" will have information about which columns contain the new
observables.The entire range of MC sweeps available inside the traj file will be
used. The second example creates an rt file as well as an averages file from the
data described in the traj file. The averages file will be called "tmp_averages".
To create histograms of the new observables, you will need to use the program
prf_his1d (see \ref prf_his1d).



\subsection usage2 Generating a sequence of PDB snapshots from a run
\verbatim
$ extract_props -op tmp_ -pdb -sf -a 4589999 -z 4600999 --every 10 n0/traj
\endverbatim

This extracts the structure of the population at cycles starting from 4589999
until 4600999 using every 10th available snapshot. The structures in pdb format
are stored in a single file (because of "-sf") called "tmp_frames.pdb".
Different snapshots in the PDB file will be given different MODEL tags.

\sa extract_snapshot, \ref settings_obs
*/

class extract_props
{
public:
    extract_props();
    ~extract_props();
    int my_init(int argc, char *argv[]);
    int parseCommand(InstructionString s);
    void show_basic_help();
    int BrowseConfs();
private:
    void auto_track_obs();
    PopulationHandler PH;
    FFHandler ffh;
    ObsHandler H;
    Etot esum;
    prf_utils::ProgArgs optn;
    unsigned long irtstart,irtend,every;
    int NTMP, log_level;
    bool NTMP_known;
    std::string mytraj,myprefix;
    prf_traj::Trajectory traj;
    prf_traj::Trajectory::iterator cb;
};
