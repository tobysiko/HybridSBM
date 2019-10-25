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

#ifndef ExVolBase_HH
#define ExVolBase_HH
#include <cmath>
#include <valarray>

namespace prf
{
//! Helper class for both ExVol and LocExVol classes.

    class ExVolBase
    {
    public:
        ExVolBase();
        virtual ~ ExVolBase();
    protected:
        int NA;
        double cut, cut2;
        double ksa;
        std::valarray < double >sigsa;
    };
}

#endif
