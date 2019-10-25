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

#ifndef SIMANNEAL_HH
#define SIMANNEAL_HH
#include "GMC.hh"

namespace prf {
    class SimAnneal : public GMC
    {
    public:
        SimAnneal();
        ~SimAnneal();
        int Setup();
        unsigned SwitchTemp();
        int parseCommand(InstructionString s);
        void print_setup();
        inline size_t n_cycles_at_T(size_t i) { return ncbeta[i]; }
    private:
        int make_default_schema();
        int read_schema();
        std::string schemafile;
        std::vector<size_t> ncbeta;
        size_t nswatT, counter;
        bool useschemafile;
    };
}

#endif // SIMANNEAL_HH
