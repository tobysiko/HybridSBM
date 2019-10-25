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

#ifndef Ran3nEngine_HH
#define Ran3nEngine_HH

#include "RandomNumberBase.hh"

namespace prf
{
    //! Ran3nEngine is the default random number generator class for profasi.
    /**
     * Ran3nEngine inherits from RandomNumberBase and implements a  method
     * given in Numerical Recipes in C. It overrides the shoot() function
     * which @returns a random number between 0 and 1.
     */

    class Ran3nEngine : public RandomNumberBase
    {
    public:
        Ran3nEngine();
        ~Ran3nEngine();
        double shoot();
        void ResetDefaultState(long seedd);
    public:
        static int inext,inextp;
        static long ma[55];
        static const long MBIG;
        static const long MSEED;
        static double FAC;
        static const long MZ;
    };
}

#endif



/* Based on method given in Numerical Recipes in C ... */
/******************************************************************/





