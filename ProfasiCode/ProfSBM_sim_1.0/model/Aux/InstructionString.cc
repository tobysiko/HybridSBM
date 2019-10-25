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

#include "InstructionString.hh"
#include "fileutils.hh"
#include "profasi_io.hh"
#include <fstream>

InstructionString::InstructionString() : mytype(COMMENT)
{
}

InstructionString::InstructionString(std::string s) : mytype(COMMENT)
{
    str(s);
}

void InstructionString::str(std::string s)
{
    s=prf_utils::trim_str(s);
    parts.clear();
    if (s.empty() or s[0]=='#') mytype=COMMENT;
    else if (s[s.size()-1]=='\\') mytype=INCOMPLETE;
    else {
        prf_utils::split(s,parts);
        if (parts[0]=="quit" or parts[0]=="exit") mytype=QUIT;
        else mytype=COMMAND;
    }
    orig=s;
}

void InstructionString::append(std::string s)
{
    if (orig[orig.size()-1]=='\\') {
        orig[orig.size()-1]=' ';
        str(orig+s);
    }
}

InstructionString InstructionString::tail()
{
    std::string rest;
    for (size_t i=1;i<parts.size();++i) rest+=(parts[i]+" ");
    return InstructionString(rest);
}

//int main()
//{
//    InstructionString s;
//    std::string line;
//
//    do {
//        if (s.type()!=INCOMPLETE) prf::cout<<"prompt> ";
//        else prf::cout<<"continuation> ";
//        getline(std::cin,line);
//        if (s.type()==INCOMPLETE) s.append(line); else s.str(line);
//        if (s.type()==COMMAND) {
//            prf::cout<<"Command : "<<s.head()<<" with "<<s.n_parts()<<" parts.\n";
//        }
//    } while (s.type()!=QUIT);
//    return 0;
//}

