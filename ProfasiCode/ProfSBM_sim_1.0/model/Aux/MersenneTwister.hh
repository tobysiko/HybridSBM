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

#ifndef MersenneTwister_HH
#define MersenneTwister_HH

#include "RandomNumberBase.hh"

const int mt_N=624;
const int mt_M=397;

namespace prf
{
    //! MersenneTwister is the default random number generator for profasi.
    /**
     * MersenneTwister inherits from RandomNumberBase. The code here is a
     * simple adaptation of the code found in GSL: The GNU Scientific Library.
     * Essentially the same working code made into a C++ class suitable for
     * use with PROFASI. One slightl change however is that we do not allow
     * the random number to be strictly zero.
     * \ingroup random_number
     */

    class MersenneTwister : public RandomNumberBase
    {
    public:
        MersenneTwister();
        ~MersenneTwister();
        double shoot();
        void ResetDefaultState(long seedd);
        void saveState(std::string sttfile);
        int recoverState(std::string sttfile);
    public:
        unsigned long mt[mt_N];
        int mti;
    };
}

#endif



/* Based on method given in Numerical Recipes in C ... */
/******************************************************************/





