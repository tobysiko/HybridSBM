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

#include "VoidEG.hh"
using std::deque;
using std::pair;

using namespace prf;
VoidEG::VoidEG() : EndGroup(VOIDEG)
{
    nd=na=0;
}

VoidEG::~VoidEG() {}

void VoidEG::Initialize() {}

void VoidEG::Reconstruct() {}

void VoidEG::Write() {}

void VoidEG::WritePDB(int &istatm, char ch_id, int aaindx, FILE * fp) {}

void VoidEG::WritePDB2(int &istatm, char ch_id, int aaindx, FILE * fp) {}

void VoidEG::Allocate() {}

void VoidEG::atomOffset(int i) {}

bool VoidEG::pep_bond_link() const
{
    return false;
}
