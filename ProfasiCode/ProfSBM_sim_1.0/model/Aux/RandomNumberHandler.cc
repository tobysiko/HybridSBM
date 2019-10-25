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

#include "RandomNumberHandler.hh"
#include "MersenneTwister.hh"
#include "Ran3nEngine.hh"
#include "profasi_io.hh"

using namespace prf;

RandomNumberHandler::RandomNumberHandler() : HandlerBase()
{
    rng=new MersenneTwister();
    explicitly_seeded=false;
    par.option("random_number_generator","rangen",
               1,"(MersenneTwister/Ran3n)");
    par.option("random_number_seed","random_number_seed",1,"(integer value)");
}

/**
\page settings_random_number Commands for the random number generator
\li \b --random_number_seed or \b -rs : Seed to be used for random number
generator initalization. If this option is not specified, a fairly random
seed is automatically generated from the UNIX PID and current time. So,
use it only if you want to retrace the same trajectory. Example: -rs 2371
\li \b random_number_generator or \b -rangen : Choose your random number
generator from Mersenne Twister and Ran3n. Example: -rangen Ran3n

*/

RandomNumberHandler::~RandomNumberHandler()
{
    if (rng!=NULL) delete rng;
}

bool RandomNumberHandler::known_generator(std::string nm)
{
    return (nm=="MersenneTwister" or nm=="Ran3n" or nm=="Ran3n");
}

int RandomNumberHandler::create_generator(std::string nm)
{
    Logger blog(10);
    if (rng!=NULL) {
        if (nm!=rng->Name()) {
            if (known_generator(nm)) {
                delete rng;
                rng=NULL;
            } else {
                blog<<"RandomNumberHandler::create_generator()> Ignoring request "
                        <<"to create unknown generator "<<nm<<"\n";
            }
        } else {
            blog<<"RandomNumberHandler::create_generator()> Already using "
                    <<nm<<"\n";
            return 1;
        }
    }
    if (nm=="Ran3n") rng = new Ran3nEngine();
    else if (nm=="MersenneTwister") rng = new MersenneTwister();
    else {
        prf::cerr<<"RandomNumberHandler::create_generator()> Should never reach "
                <<"this line!\n";
        return 0;
    }
    return 1;
}

void RandomNumberHandler::set_seed(int sd)
{
    rng->setSeed(sd);
    explicitly_seeded=true;
}

void RandomNumberHandler::auto_seed(int rnk)
{
    rng->RandomizeState(rnk);
}

int RandomNumberHandler::parseCommand(InstructionString s)
{
    if (s.head()=="random_number_generator") {
        create_generator(s.tail().str());
    } else if (s.head()=="random_number_seed") {
        set_seed(atoi(s.tail().str().c_str()));
    }
    return 1;
}
