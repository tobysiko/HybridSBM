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

#ifndef VoidEG_HH
#define VoidEG_HH
#include "EndGroup.hh"

namespace prf
{

    class VoidEG : public EndGroup
    {
    public:
        VoidEG();
        ~VoidEG();
        void Allocate();
        void atomOffset(int i);
        bool pep_bond_link() const;
        void Initialize();
        void Reconstruct();
        void Write();
        void WritePDB(int &istatm, char chid, int aaind,
                      FILE * fp);
        void WritePDB2(int &istatm, char chid, int aaind,
                       FILE * fp);
    };
}

#endif
