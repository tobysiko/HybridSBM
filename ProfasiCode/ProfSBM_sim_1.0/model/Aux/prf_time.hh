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

#ifndef PRF_TIME_HH
#define PRF_TIME_HH
#include <ctime>
#include <sys/time.h>
#include <string>

namespace prf_utils
{
    //! Abstraction for some system time function
    class prf_time
    {
    public:
        prf_time();
        ~prf_time();
        prf_time(const prf_time &);
        prf_time & operator=(const prf_time &);
        void update();
        //! Time interval since some given time in seconds
        double operator-(const prf_time &);
        std::string to_UTC();
        std::string now();
        std::string stamp();
    private:
        // The internal representation of time is not important for
        // PROFASI. All that matters is a way to tell what time it is
        // now, and how long it has been since X. There are many
        // system functions to find the current time. To avoid
        // getting the mess of time_t, timeval etc into PROFASI,
        // this class provides a useful interface, and hides the
        // details.
        timeval t;
    };
}

#endif
