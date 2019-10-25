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

#ifndef BasicMCRun_HH
#define BasicMCRun_HH
#include <Elements/PopulationHandler.hh>
#include <Energy/FFHandler.hh>
#include <Observables/ObsHandler.hh>
#include <Algorithms/MC.hh>
#include <Aux/RandomNumberHandler.hh>
#include <Aux/Trajectory.hh>

using namespace prf;
using std::map;
using std::string;

/**
  \page bmcrun_progref BasicMCRun program reference
The program BasicMCRun performs canonical Monte Carlo simulations on single
or multiple chain systems at a constant temperature using the ProFASi libraries.
It also forms the base for  programs like simulated annealing, simulated
tempering, parallel tempering.

Example use:
\verbatim
BasicMCRun --add_chain 1 "*GEWTY DDATKT FTVTE*" -ncyc 100000 -nrt 1000 -T "300 Kelvin"
\endverbatim

In the above example, BasicMCRun was asked to perform a short simulation of 100000
MC sweeps for a 16 residue peptide (residues 41--56 of 1GB1.pdb), at 300 Kelvin.
The program was asked to write to the run-time history file "rt" every 1000 sweeps.

The progam accepts a large number of options which can be given as command line
arguments or in a settings file. This file is called "settings.cnf" by default,
although this name can be changed using a command line option "--settings_file
yourfilename". All command line options in their long form can be used as
settings file instructions. For instance, the above run could have been
done with the following settings.cnf file:
\verbatim
add_chain 1 * GEWTY DDATKT FTVTE *
num_cycles 100000
rt_write_freq 1000
temperature 300 Kelvin
\endverbatim
With a file called "settings.cnf" with the above content, one can start BasicMCRun
without any command line arguments with exactly the same effect as in the first
example above. Notice that the quote marks in the command line example above
have been dropped in the settings file. Also, those familiar with older versions
of ProFASi will recognize that the sequence is no longer enclosed in
angle brackets &lt;&gt;. It does not hurt to write the sequence in the
old way, i.e., "add_chain 1 <* GEWTY DDATKT FTVTE *>", but it also does not
do any good. In older versions of ProFASi, sequence input continued until the
termination of the angle brackets, so that the sequence could be written
in many lines. Now, sequence reading continues till the end of the line.
If the sequence is long, and you want to write it in many line, you must
use a continuation character "\", as shown below.
\verbatim
add_chain 1 * GEWTY \
              DDATKT \
              FTVTE *
num_cycles 100000
rt_write_freq 1000
temperature 300 Kelvin
\endverbatim
Notice also that in the settings file, we use the long forms of the commands:
"temperature" instead of "T". On the command line, the temperature can be
specified using --temperature "300 Kelvin" or -T "300 Kelvin". But only the
long form is recognized in the settings file.

The complete list of options accepted by the program is very long. But they
can be grouped into a few categories. BasicMCRun accepts
\li \ref settings_population
\li \ref settings_random_number
\li \ref settings_mc
\li \ref settings_obs

In addition, the following miscellaneous instructions control various aspects
of a run with BasicMCRun.

\li \b --num_cycles or \b -ncyc: Total number of Monte Carlo sweeps (not steps) to be performed.
Example: -ncyc 10000000. Aliases (only valid in settings file): ncycles,nsweeps,mccyc
\li \b --num_therm_cycles or \b -ntherm : The number of cycles at the start of the run that are
discarded when making histograms or calculating averages. Example: -ntherm 30000
\li \b rt_write_freq or \b -nrt :  Sets the frequency at which system properties are recorded in
the run-time history file "rt". It is measured in the unit of cycles. A value of 1 means one Monte
Carlo cycle and not one Monte Carlo step. Example: --rt_write_freq 1000
\li \b --conf_write_freq or \b -iconf: Frequency of saving system configuration.
Example: -iconf 1000
\li \b --avg_write_freq or \b -iavg : Frequency of calculating averages.
Example: -iavg 10000
\li \b --num_progress_reports or \b -nreports : The number of times the program tries to
estimate how much longer it has to run to finish the job. If the estimated time-requirement
exceeds the remainder of the allocated time, it estimates what MC cycle it will reach when the
time is up. Example: --num_progress_reports 10
\li \b --log_level or \b -ll : Control the verbosity level of the log messages.
Example: --log_level 10
\li \b --secondary_settings or \b -secondary_settings : Whether there are node specific
settings specified in directories n0,n1 ... The secondary settings files have to be called
"settings.cnf" and must be found in the node specific directories n0, n1 etc. The secondary
settings mechanism is still supported in ProFASi v. 1.5, but it is deprecated. To issue
rank specific instructions, use the command "for_rank ..." in the top level
settings file (See description in \ref secondary_settings).
Example: --secondary_settings or --no-secondary_settings
In a settings file, you would write this command as "secondary_settings on".
\li \b --stdout_redirect or \b -stdout_redirect: If set, each node writes standard output
messages to a separate file n0/stdout, n1/stdout etc.
Example: --stdout_redirect or --no-stdout_redirect
\li \b --stderr_redirect or \b -stderr_redirect: If set, each node writes standard error
messages to a separate file n0/stderr, n1/stderr etc.
Example: --stderr_redirect or --no-stderr_redirect


    \sa \ref settings_random_number, \ref settings_population, \ref settings_mc,
    \ref settings_obs

*/

//! Constant temperature Monte Carlo simulations
/**
This is an application program written in the form of a class. It has a Run
function, and the "main" function simply needs to transfer control to this Run
function. In ProFASi, there is no such thing as the "main profasi  executable".
Application programs should be written in the form of an interface class
like this one and linked to the ProFASi libraries.
\sa \ref bmcrun_progref
  */
class BasicMCRun
{
public:
    BasicMCRun();
    virtual ~BasicMCRun();
    // Return value of 0 means initialisation failed. 1 means success.
    // Return value of -1 means the program should quit without doing
    // anything, but yet, without printing more error messages. It is
    // then assumed that some function call during Init has already
    // printed the error message.
    int init(int argc,char *argv[]);
    /*
      * In ProFASi version 1.5, command parsing happens in many stages.
      * Each class, like BasicMCRun, only parses the "lightest" commands
      * immediately, like for instance, commands asking to change the
      * value of some parameter. More complex commands such as those
      * changing sequence, updates, or even the force field are collected,
      * and executed when appropriate.
      */
    //! Execute commands collected with getCommands
    int parseCommands(std::list<InstructionString> &cmds);

    virtual int parseCommand(InstructionString s);
    void disable(std::string anoption);
    int Run();
public:
    inline void set_rank_nruns(int i1,int i2) {myrank=i1;nruns=i2;}
protected:
    void init_streams();
    int init_dir();
    virtual void init_filenames();
    virtual int init_MC();
    virtual int get_span(prf_traj::Trajectory &traj,
                         unsigned long &i1, unsigned long &i2);
    virtual void write_MC_averages();
    void print_MC_setup();
    virtual void writeTemperatures();
    virtual void run_relaxation_cycles();
    virtual void help(std::string qury);
    virtual void show_help();
    virtual void show_basic_help();
    //The following are relevant mostly for derived classes
    //indx=temperature index,icycl=cycle number. Called when a run is to be
    //reset to a certain point in a configuration file
    virtual void set_index(int indx, unsigned long icycl);
    virtual int update_T();
    // In BasicMCRun should_suspend(bool) trivially returns the input value
    // In ParTempRun it checks across all runs if any of them want to suspend
    // and goes by the gloomiest prediction among all replicas
    virtual bool should_suspend(bool myrec);
    void auto_track_obs();
    void printMCconfig();
    //! Write program state to PROFASI binary configuration file
    /**
      * Writes information about the current state of the program to a
      * compressed binary configuration file. Every time this function
      * is called, the following information is written:
      * (i) MC time (ii) Temperature index (This is always
      * 0 for BasicMCRun, but contains useful information for programs
      * like ParTempRun. Writing a redundant 0 in BasicMCRun keeps the
      * format of the data in the binary files fixed across all simulation
      * programs in PROFASI) (iii) Total energy (iv) State of the
      * random number generator, normally, as a pair of long and unsigned
      * long variables. The first is the original seed, the second is the
      * number of calls to the generator up until the time of the
      * configuration save. Using these two, the random number generator
      * can be restored (v) All degrees of freedom of the population.
      *
      * The data is simply dumped without any formatting. So, the contents
      * of the state files (normally called "conf...", found inside the
      * run directories) are not (normal-) human readable. Since this way
      * of writing data makes use of the internal structure of the data,
      * it achieves optimal compression. But without knowledge of the
      * layout of the bytes in the files, this data would be totally
      * unreadable. Such information was formerly contained in a separate
      * file called "conf.info" also generated by the simulation programs.
      * The same byte for byte description of the contents of each block of
      * data in the conf files is now contained as a human readable "header"
      * section inside the "conf..." files.
      *
      * PROFASI 1.5 makes the interaction with the binary "conf..." files
      * more indirect, through the PROFASI trajectory description file,
      * called n0/traj, n1/traj etc. This is an important topic and is
      * described in detail in \ref traj_files .
      */
    virtual void writeConf();
    int close_traj_file();
    //! Restore the state of the run to a given point of a previous run
    /**
      * This restores the population, random number generator, temperature
      * and energy classes to a specified point in a previously generated
      * trajectory. MC updates should then follow the path of the older run
      * unless external factors change this. The parameter passed is the
      * MC cycle where the program should be restored. If a large value
      * is passed (or the default argument used), the run is restored to
      * the last saved configuration in the older trajectory. The name of the
      * trajectory file used is fixed: it's called "traj", and it is located
      * in the output directory of the run. Since the "traj" files are only
      * envelopes containing information about data (but no actual data),
      * additional files are needed for the recovery to work. These are the
      * conf.data___ files found in the output directories. Without the
      * conf.data___ files, the traj files can not be used to restore the
      * program.
      */
    int recover(unsigned rcyc=(unsigned)-1);
    //! Open a new conf.data___ file for a new simulation segment
    /**
      * Open a data file to contain binary configuration snapshots for the
      * run. Every time a run is restarted, a new conf.data___ file is created
      * with a name based on the time of the first write event. The name is
      * registered inside the "traj" file for the run, so that the traj file
      * contains information about the segments in the appropriate order.
      */
    FILE *open_segment(std::string prfx);

    int curTindex,NTMP;
    ProgArgs optn;
    unsigned int NRT,IAVG,ICONF,NTHERM,INORM,NTESTI,rem_progress_writes;
    int snapshot_format,log_level;
    int ncyc_per_T_updt;
    unsigned long MCCYC, RMCCYC, icyc, icyc0, nblk;
    std::string locdir,pdbfil,minpdb,statfile,speedfile;
    std::string trajfile,rnboosterfile,datafl;
    // Program switches
    std::map<string,bool> swtch,usrspc;
    // Run specific variables for help with parallelization
    int myrank,nruns;
    string mydir, progheader,progname;
    // Miscellaneous
    Etot esum; //total energy
    Ego ego; // native contact energy
    std::list<string> postponed;
    MC canon;
    MC *mc;
    double Tcur, Trelax, available_run_time,en;
    //! A population handler
    PopulationHandler PH;
    //! An observable handler
    ObsHandler H;
    //! A handler object for random numbers
    RandomNumberHandler ranh;
    //! A force field handler
    FFHandler ffh;
    
    bool debug; // set me in the constructor!
};

#endif
