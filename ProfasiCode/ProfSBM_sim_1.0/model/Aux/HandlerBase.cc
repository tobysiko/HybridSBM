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

#include "HandlerBase.hh"
#include "profasi_io.hh"
HandlerBase::HandlerBase() {}

HandlerBase::~HandlerBase() {}

void HandlerBase::parseCommands(std::list<InstructionString> &cmds, int argc, char *argv[])
{
    par.analyze(argc, argv);
    par.add_cmd_line_instructions();
    std::list<InstructionString> myextras=par.get_options();
    for (std::list<InstructionString>::iterator it=cmds.begin();it!=cmds.end();++it) {
        parseCommand(*it);
    }

    for (std::list<InstructionString>::iterator it=myextras.begin();it!=myextras.end();++it) {
        parseCommand(*it);
    }
}

int HandlerBase::parseCommand(InstructionString s) {
    return 1;
}

void HandlerBase::show_help() {
    par.write_available();
}

bool HandlerBase::args_check(InstructionString s, size_t expected)
{
    if (s.n_parts()<(expected+1)) {
        prf::cerr<<"Incorrect syntax for "<<s.head()
                <<": expected "<<expected<<" arguments. \n";
        return false;
    } else return true;
}
