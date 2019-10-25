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

#ifndef INSTRUCTIONSTRING_HH
#define INSTRUCTIONSTRING_HH

#include <string>
#include <deque>
#include "fileutils.hh"
#include <fstream>

typedef enum {QUIT, COMMENT, INCOMPLETE, COMMAND} InstructionType;
class InstructionString
{
public:
    InstructionString();
    InstructionString(const std::string s);
    void str(std::string lin);
    void append(std::string lin);
    inline std::string str() const { return orig; }
    inline InstructionType type() const { return mytype; }
    inline size_t n_parts() const { return parts.size(); }
    inline std::string part(size_t i) const { return parts[i]; }
    inline std::string head() const { return parts[0]; }
    InstructionString tail();

private:
    InstructionType mytype;
    std::deque<std::string> parts;
    std::string orig;
};

template <class T>
int get_commands(std::string filename, T &cmds)
{
    if (!prf_utils::TestFile_r(filename)) return 0;
    std::ifstream fin(filename.c_str());
    return get_commands(fin,cmds);
}

template <class T>
int get_commands(std::istream &stin, T &cmds)
{
    cmds.clear();
    InstructionString s;
    std::string line;
    while (getline(stin,line)) {
        if (s.type()==INCOMPLETE) s.append(line); else s.str(line);
        if (s.type()==COMMAND) cmds.push_back(s);
    }
    return 1;
}

#endif // INSTRUCTIONSTRING_HH
