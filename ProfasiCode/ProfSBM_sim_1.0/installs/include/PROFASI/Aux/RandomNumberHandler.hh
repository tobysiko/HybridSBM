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

#ifndef RANDOMNUMBERHANDLER_HH
#define RANDOMNUMBERHANDLER_HH
#include "HandlerBase.hh"
#include "RandomNumberBase.hh"

namespace prf {
    class RandomNumberHandler : public HandlerBase
    {
    public:
        RandomNumberHandler();
        ~RandomNumberHandler();
        inline RandomNumberBase *generator() { return rng; }
        bool known_generator(std::string nm);
        int create_generator(std::string nm);
        int parseCommand(InstructionString s);
        void auto_seed(int rnk);
        void set_seed(int sd);
    private:
        RandomNumberBase *rng;
        bool explicitly_seeded;
    };
}

#endif // RANDOMNUMBERHANDLER_HH
