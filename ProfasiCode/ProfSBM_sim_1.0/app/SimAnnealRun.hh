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

#ifndef SimAnnealRun_HH
#define SimAnnealRun_HH
#include <Algorithms/SimAnneal.hh>
#include "BasicMCRun.hh"

using namespace prf;

/**
  \page sa_progref SimAnnealRun program reference
  \section sa_desc Brief description of simulated annealing
In simulated annealing, canonical Monte Carlo simulations are carried out with
a progressively decreasing temperature. At the higher temepratures the MC
simulation samples the energy landscape with a coarse resolution, and as the
temperature goes down, it starts to distinguish between finer and finer details.
If the cooling is slow enough and the process of high temperature equilibration
and the subsequent cooling is repeated many times, the global minimum in the
energy landscape is very likely to be found.

The program SimAnnealRun performs simulated annealing using ProFASi libraries.
The actual temperature at each temperature step, as well as the number of Monte
Carlo sweeps performed at that temperature can be controlled through a "schema".
The default is a geometric progression of temperatures (constant ratio between
adjacent temperatures), with a fixed number (500) of Monte Carlo sweeps at each
temperature. But any alternative schema can be provided using a schema file,
specified by the "schema_file" command in the settings file or on the command
line as described below in \ref sa_schemafile.

SimTempRun starts an annealing cycle by initialising the population in a
specified manner. For instance, if option init_config was set to "stretched",
the chains are put to a stretched conformation. If it is "random" (default), a
random conformation is used. If it is "file://someconfig.tcnf", the
population is initialised to the conformation stored in the file
"someconfig.tcnf". It then performs a specified number of Monte Carlo sweeps
at the highest temperature. Then the temperature is reduced to the next lower
temperature and the required number of sweeps at that temperature are carried
out. When the required number of Monte Carlo sweeps at all temperatures are
finished, a new annealing cycle is started at the highest temperature.

An annealing cycle is started only if there is a sufficient number of Monte
Carlo cycles remaining. For instance, if there are 20 temperatures, each with
1000 sweeps, and a total number of 50000 sweeps are requested, the run will stop
when 2 annealing cycles, i.e., 40000 sweeps are finished.

<div style="color:red">It is recommended that you use the switch
"--no-histograms" to suppress generation of histograms if you decide
to run simulated annealing with a large number of temperatures and only a few
sweeps at each temperature. </div>

\section sa_usage Usage
Most options available on the command line in the BasicMCRun program work with
SimAnnealRun. To see how to set up the population, force field etc, please
see \ref bmcrun_progref. The following are some extra options for simulated
annealing.

\li \b --schema_file or \b -sch : Name of the file containing a schema for
lowering the temperature during simulated annealing.
\li \b --min_temperature or \b -tmin : Minimum temperature in Kelvin
\li \b --max_temperature or \b -tmax : Maximum temperature in Kelvin
\li \b --num_temperatures or \b -ntmp : Number of temperatures
\li \b --cycles_at_T or \b -ncT : Number of cycles at each temperature

Note that if a schema file is specified with the "--schema_file" option, minimum
and maximum temperatures, number of temperatures and number of cycles at each
separate temperature is read from the schema file. The other options above are
then ignored.
\sa \ref bmcrun_progref

\section sa_schemafile Simulated Annealing schema file
The schema file should contain 2 columns: temperature and the number of sweeps
at that temperature. The temperature values are assumed to be in the PROFASI
internal units by default. To enter temperatures in Kelvin, you should have a
comment line at the top with the entry "#temperature Kelvin". To enter inverse
temperatures, that line should be "#temperature inverted". In order to give the
user complete freedom to devise any temperature "lowering" scheme, the
SimAnnealRun program simply follows the temperatures in the order they appear
in the schema file, without any attempt to sort them. So, the first entry should
normally be the highest temperature with successive entries at lower and lower
temperatures. You can arbitrarily choose how many sweeps to perform at each
temperature. Here is an example schema file.

\verbatim
#temperature Kelvin
400 100
350 120
330 500
320 1000
315 3000
310 1000
300 1000
290 1000
280 1000
270 1000
\endverbatim

\sa \ref bmcrun_progref

*/

//! Simulated annealing
/**
  This class implements simulated annealing using ProFASi libraries. It is
  based on the class BasicMCRun. Most properties are directly inheritted.
  Only the temperature update method is implemented here as a cooling
  step. There is also a function to reset the system to the initial state
  at the start of each simulated annealing cycle. To see how to use the
  program, see \ref sa_progref
  \sa \ref sa_progref
  */
class SimAnnealRun : public BasicMCRun
{
public:
    SimAnnealRun();
    ~SimAnnealRun();
    //! Parse commands specific to simulated annealing
    int parseCommand(InstructionString s);
    int init_MC();
    void show_basic_help();
    int update_T();
    int re_init();
    void writeTemperatures();
private:
    SimAnneal siman;
};

#endif
