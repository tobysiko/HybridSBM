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

#ifndef HANDLERBASE_HH
#define HANDLERBASE_HH
#include "ProgUtils.hh"

class HandlerBase
{
public:
    HandlerBase();
    virtual ~HandlerBase();
    //! Parse commands passed as a deque of InstructionString objects
    virtual void parseCommands(std::list<InstructionString> &cmds, int argc, char *argv[]);
    virtual int parseCommand(InstructionString s);
    virtual void show_help();
protected:
    bool args_check(InstructionString s, size_t expected);
    prf_utils::ProgArgs par;
};

#endif // HANDLERBASE_HH
