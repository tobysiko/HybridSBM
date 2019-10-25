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

#include "prf_time.hh"
#include "profasi_io.hh"
#include <unistd.h>
#include <sstream>

namespace prf_utils
{
    prf_time::prf_time()
    {
        update();
    }

    prf_time::~prf_time() {}

    prf_time::prf_time(const prf_time &gt) : t(gt.t) {}

    prf_time & prf_time::operator=(const prf_time &gt)
    {
        if (this!=&gt) {
            t=gt.t;
        }

        return *this;
    }

    void prf_time::update()
    {
        gettimeofday(&t,NULL);
    }

    double prf_time::operator-(const prf_time &gt)
    {
        return difftime(t.tv_sec,gt.t.tv_sec)+
               (1e-6)*(t.tv_usec-gt.t.tv_usec);
    }

    std::string prf_time::to_UTC()
    {
        char strtime[200];

        struct tm tmptime;
        gmtime_r(&t.tv_sec, &tmptime);
        strftime(strtime,sizeof(strtime),"%Y-%h-%d-%H:%M:%S",
                 &tmptime);
        return std::string(strtime);
    }

    std::string prf_time::now()
    {
        update();
        std::ostringstream ost;
        ost<<t.tv_sec;
        return ost.str();
    }

    std::string prf_time::stamp()
    {
        update();
        char outstr[200];
        struct tm tmptime;
        gmtime_r(&t.tv_sec,&tmptime);
        strftime(outstr,sizeof(outstr),"%Y-%h-%d-%H%M%S",&tmptime);
        return std::string(outstr);
    }

}
