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

#ifndef MINIMIZE_HH
#define MINIMIZE_HH
#include <Elements/PopulationHandler.hh>
#include <Energy/FFHandler.hh>
#include <Observables/ObsHandler.hh>
#include <Observables/ObsExpression.hh>
#include <Algorithms/Minimizer.hh>

//! Minimize an expression over degrees of freedom of a protein molecule
/**
  The program minimize uses conjugate gradient minimization to minimize a given
  observable or a linear combination of observables near a given structure.
  <i>Typically this program is used with both settings file and command line
  options together.</i> For instance, if we have the following settings file:
  \verbatim
  set_population testprot.pdb::A
  new_obs ProteinRMSD rmsd using +HV ; struc1 testprot.pdb::A ; struc2 $::A
  new_obs Rg rg
  \endverbatim

  One can then minimize an arbitrary linear combination of the defined
  observables, like this:
  \verbatim
  minimize "Etot+25*rg+10*rmsd"
  \endverbatim
  \e Etot in the above example is the total energy, which is an implicitly
  defined observable, even if it is not created with a "new_obs" command.
  Similarly, the component terms of the energy function can be used without
  their explicit declaration in the settings file.

  Minimization can be restricted to only the side chain degrees of freedom with
  the option \e --using.
  \verbatim
  minimize "Etot+25*rg+10*rmsd" --using "sc" "*"
  \endverbatim
  The \e --using option takes 2 parameters. It is <i> --using what where </i>.
  The first argument could be "bb", "sc" or "all", for backbone, side chain and
  all degrees of freedom. The second argument is a selection inside the
  population where the degrees of freedom will be varied. For instance,
  \verbatim
  minimize "rg" --using "all" "A/15/56"
  \endverbatim
  will minimize the observable \e rg using all degrees of freedom in chain A
  for group indices (start from 0, include capping groups) 15 through 56.

  The output of the minimize program is stored in a file called \e output.xml
  in \ref xmlstruc .
  \sa \ref regul , \ref xmlstruc, \ref toc

  */
class minimize
{
public:
    minimize();
    ~minimize();
    int init(int argc, char *argv[]);
    int run();
    void help(std::string tpc);
    void show_basic_help();
    void parseCommands(std::list<InstructionString> &cmds);
    void parseCommand(InstructionString &s);
    void auto_track_obs();

private:
    PopulationHandler PH;
    FFHandler ffh;
    Etot esum;
    ObsHandler H;
    Minimizer cg;
    ObsExpression m;
    ProgArgs optn;
    std::string outputfile, oformat, expression,sel,typ;
    std::vector<DOF_Info> alldofs;
};

#endif // MINIMIZE_HH
