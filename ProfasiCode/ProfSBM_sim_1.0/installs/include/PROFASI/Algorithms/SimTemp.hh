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

#ifndef SIMTEMP_HH
#define SIMTEMP_HH
#include "GMC.hh"

namespace prf {
    class SimTemp : public GMC
    {
    public:
        SimTemp();
        ~SimTemp();
        int Setup();
        unsigned SwitchTemp();
        void read_g_pars(std::string gfl);

        //! Write recommendation for a new set of g parameters to file
        void writeNewGPars(std::string flnm);
        int parseCommand(InstructionString s);
        void print_setup();
    private:
        std::string gpfile;
        std::vector<double> g;
        size_t nswatmpt;
    };
}
#endif // SIMTEMP_HH
