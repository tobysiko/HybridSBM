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

#ifndef STAN_HH
#define STAN_HH

#include "ObsHandler.hh"
#include <Aux/InstructionString.hh>

using namespace prf_app;

class StAn
{
public:
    StAn();
    ~StAn();
    int load(std::string filename);
    void pre_init();
    void re_init_obs();
    void save(std::string filename);
    Matrix<double> profile(std::string exprs,std::string rng);
    void interactive_session();
    void process_command(InstructionString s);
//    void process_commandline();
    void show(InstructionString s);
    void parseCommands(std::string flnm);
private:
    bool composition_changed();
    Population p;
//    std::vector<prf::Energy *> eterm;
    ObsHandler H;
    std::vector<std::vector<prf::OneLetterCode> > thecomposition;
};

#endif // STAN_HH
