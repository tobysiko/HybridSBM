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

#ifndef WLRun_HH
#define WLRun_HH
#include "BasicMCRun.hh"
#include <Algorithms/WangLandau.hh>
#include <Aux/ProgUtils.hh>

using namespace prf;

//! Monte Carlo simulations with the Wang Landau method
/**
  * This class inherits from the BasicMCRun class and provides the
  * necessary virtual function overrides, so that the simulations actually
  * perform Wang Landau iterations instead of Metropolis-Hastings Monte Carlo.
  * It uses large parts of the code in BasicMCRun, including its Run and
  * Init functions. The MC pointer is assigned to a prf::WangLandau object,
  * which is where the details of the method are implemented.
  */
class WLRun : public BasicMCRun
{
public:
    WLRun();
    ~WLRun();
    // Parse commands which must be processed in this class
    int parseCommand(InstructionString s);
    void show_basic_help();
    void writeConf();
protected:
    int init_MC();
    void write_MC_averages();
    void run_relaxation_cycles();
    int update_T();
private:
    WangLandau wl;
    int maxfupdates, nwlbins;
};

/**
  \page wlrun_progref Wang-Landau method with WLRun
The program WLRun performs Wang-Landau iterations. The basic algorithm is as
described in "Fugao Wang and D.P. Landau, Physical Review Letters,
86 (10) 2050", with the following modifications:

\li If the energy exceeds the given upper limit, an update is by default
rejected. But if the flag --use_htmc is given, the update is accepted or
rejected as if it was a Metropolis MC update with a high temperature. New
Metropolis MC updates are done until the energy comes back into the
specified range again.
\li For protein systems, it is not possible to know ahead of time what
a good lower bound for energy should be. If the lower limit is set below
what is possible for the given system, the simulation will never converge.
Therefore, the simulations should be started with a conservative estimate
of the lower bound. It is ok if it is too high. If the simulation finds a
state with energy lower than the lower bound, it is not rejected. It is
kept, and the lower bound is adjusted by adding new bins until the range
includes the newly found state. If the range is changed, the density of
states must now be determined over a wider range. The simulation continues
with the existing estimate of the density of states for the pre-existing
bins.

As of now, the Wang Landau algorithm is less well tested in ProFASi than
the other simulation algorithms. Please carefully check the code and
test it thoroughly before using it for anything serious. Of course,
if you find a bug or improve its performance, we would like to know.

\section options Options
This implementation uses the program BasicMCRun as a base class. Most
options available for BasicMCRun can be used for WLRun. In particular
setting up the protein population, updates etc. proceeds exactly as
in BasicMCRun. In addition to those options, this program accepts
\li \subpage wlopts
\section Example
\verbatim
$ WLRun --energy_range 20 70 --n_bins 50 -htmc --n_stages 10 \
  --add_chain 1 "*GEWTYDDATKTFTVTE*" -ncyc 1000000 -lcyc -1 -time 9:00:00
\endverbatim
This is how the above configures a Wang-Landau run:
\li The energy histogram is initially between 10 and 70 profasi units
\li The histogram initially has 50 bins. It will grow when lower energy states
are found.
\li It will use the high temperature MC method to put back a trajectory which
exits the given range on the higher side.
\li It will try to perform 10 rounds of refinement of the Wang-Landau update
factor f: f-->sqrt(f).
\li The population will contain 1 chain of sequence "*GEWTYDDATKTFTVTE*". The
initial and final asterix ("*") symbols are toggle markers in ProFASi to switch
back and forth between 3 letter codes and one letter codes for residues. The
same sequence could be written without the asterixes as "GLY GLU TRP THR TYR
ASP ASP ALA THR LYS THR PHE THR VAL THR GLU".
\li The run will have 1000000 Monte Carlo sweeps, unless it converges before
that.
\li \b -lcyc \b-1 : Each sweep will consist of as many MC updates as there are
degrees of freedom in the system. You could write -lcyc 100 to force 100 updates
per sweep. The value -1 is a convenient shorthand for the number of degrees of
freedom in the current system.
\li The run will suspend automatically if it is not finished within 9 hours.

*/

#endif
