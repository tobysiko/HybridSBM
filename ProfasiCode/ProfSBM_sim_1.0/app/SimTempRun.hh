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

#ifndef SimTempRun_HH
#define SimTempRun_HH
#include <Algorithms/SimTemp.hh>
#include "BasicMCRun.hh"

using namespace prf;

/**
  \page simtemp_progref SimTempRun program reference

This program performs simulated tempering runs on single chain and multi chain
systems.

\section simtemp Brief description of simulated tempering
Simulated tempering is a generalised ensemble Monte Carlo method in which
the temperature of the Monte Carlo is regarded as as dynamic variable.
This means, during a simulation the temperature can
vary in a given discrete set of values  (t<sub>0</sub>,t<sub>1</sub>,
t<sub>2</sub>, ...t<sub>n-1</sub>), with t<sub>0</sub> &gt; t<sub>1</sub>
&gt; t<sub>2</sub> &gt; ... t<sub>n-1</sub> . The inverse temperatures
\f$\beta_i\f$ satisfy \f$\beta_0 < \beta_1 < ...\f$ By default, this discrete
set of temperature values is chosen to be a geometric series: t<sub>i</sub>
= t <sub>0</sub> r<sup>i</sup> . After a predetermined number of elementary MC
updates (typically after every sweep, which consists of N<sub>dof</sub>
elementary MC updates in ProFASi v. 1.5 by default), an update of the
temperature variable is attempted. The system temperature can change
to one of its nearest neighbouring temperatures in the temperature set, and
this is a random event with probability depending on total energy E, the
temperatures and a set of simulated tempering parameters
(g<sub>0</sub>, g<sub>1</sub>,...,g<sub>n-1</sub>) corresponding to the
temperatures: \f$\exp(-(\Delta \beta E + \Delta g))\f$. \f$\Delta \beta \f$
and \f$\Delta g\f$ are the changes of (inverse) temperature and g.

The g-parameters are unknown properties of the system, which must be determined
through the simulations. This is done using an iterative procedure as will be
explained now. We run a simulated tempering simulation with a
certain set of g-parameters, and observe how often different temperatures in
the temperature set are visited. We calculate a new set of g-parameters,
given by g'<sub>i</sub> = g<sub>i</sub> + log(p<sub>i</sub>), where
p<sub>i</sub> are the probabilities of occurence of different temperatures.
The g' values are then used as the g-parameters for a new simulation. This
process is repeated until the calculated changes to the g-parameters are
negligible. Tuturial \ref tutorial_3 gives a walk through of this procedure
with ProFASi.

\section simptemp_proguse Usage
Except for the multiple temperatures and temperature updates, simulated
tempering is just like a canonical Monte Carlo. So, this program is also
just like BasicMCRun in its usage. All options that can be passed to that
program can also be used with this. To see how to set up the protein population,
the force field, random numbers etc, look at the \ref bmcrun_progref. In
addition there are a few options unique to simulated tempering, which will be
described here.

\li \b --ntmp or \b -ntmp 8 : Number of temperatures
\li \b --tmin or \b -tmin : "274 K" or -tmin "274 Kelvin" Lowest temparature.
If "Kelvin" or "K" is not written, the internal PROFASI units will be used.
\li \b --tmax or \b --tmax : 374 K Highest temperature
\li \b --g_parameter_file \b -g_parameter_file : gpars.in. Input file with
simulated tempering g-parameters. For a given system, the parameters in this
file are successively refined after each run.
*/

//! Simulated tempering with profasi
/**
The class definition here looks a bit thin. It does not even have an explicit
Run function. But it inherits a large number of properties from the class
BasicMCRun. Most of the steps required during a simulated tempering run are
also required in a simple cannonical MC run. There are only a few extra steps.
This class implements only those differences relative to BasicMCRun. In
simulated tempering, the temperature is a dynamic variable, and changes during
the simulation. The algorithm defines a specific way in which the temperature
is updated, which is implemented in the update_T() function. The Run function
of the BasicMCRun class calls an update_T function at a suitable point. That
has no meaning for a cannonical MC run, but gives a way to implement other
algorithms without repeating too much code. The parseCommands function here
only processes commands that are specific to simulated annealing.
 */

class SimTempRun : public BasicMCRun
{
public:
    SimTempRun();
    ~SimTempRun();
    void show_basic_help();
protected:
    int init_MC();
    void init_filenames();
    int update_T();
    void write_MC_averages();
    void writeTemperatures();
private:
    SimTemp stmp;
    std::string ogpfile, gpstfil;
};

#endif
