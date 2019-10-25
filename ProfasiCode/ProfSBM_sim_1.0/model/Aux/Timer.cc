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

#include "Timer.hh"
#include "profasi_io.hh"
#include <unistd.h>

namespace prf_utils
{
    Timer::Timer() : thegoal(10000000), thestart(0), wherenow(0),
            projected_pos(0), themaxtime(1e20), treq(0),tavail(1e20),
            npt(10), limited(false), mission_abandonned(false),
            advice(go_on), posv("Position") {}

    Timer::~Timer() {}

    void Timer::Start(unsigned long startpos)
    {
        milestones.resize(npt+1,0);
        mszone=(thegoal-thestart)/(double) npt;

        for (int i=0;i<npt+1;++i) milestones[i]=(unsigned long)((i+1)*mszone);

        mszone=0.5*mszone;

        ipt=0;

        time_strt.update();

        tnow=time_strt;

        projected_pos=wherenow=thestart=startpos;

        treq=0;

        rec_interval=tavail=themaxtime;

        advice=go_on;

        expecting_finish=report_asap=true;

        mission_abandonned=false;

        prf::cout<<"Timer started with "<<posv<<" "
        <<wherenow<<" at UTC = "<<tnow.to_UTC()<<"\n";
    }

    void Timer::Stop(unsigned long endpos)
    {
        wherenow=endpos;
        time_end.update();
        tnow=time_end;

        if (goal_reached()) {
            prf::cout<<"Job finished.\n";
            PrintInterval();
        } else if (mission_abandonned) {
            prf::cout<<"The final "<<posv<<" reported to the timer is "
                    <<wherenow<<". The original target was "<<thegoal
                    <<", but the timer has been told that we are not "
                    <<"going any further right now.\n";
            PrintInterval();
        } else {
            prf::cout<<"Job unfinished. Final "<<posv<<" reported "
            <<" to the timer is "<<wherenow<<". The goal was "
            <<thegoal<<". Starting from "<<wherenow
            <<", it will take "<<Seconds(treq)
            <<" more to finish the requested task.\n";
        }
    }

    double Timer::Measurement()
    {
        time_tmp.update();
        return (time_tmp-time_strt);
    }

    void Timer::PrintInterval()
    {
        double tsec=(time_end-time_strt);
        prf::cout<<"Time Interval measured by the timer is : "
        <<Seconds(tsec)<<"\n";
    }

    void Timer::Reset()
    {
        time_tmp.update();
        time_end=time_strt=time_tmp;
    }

    void Timer::record(unsigned long progr)
    {
        unsigned long wherethen=wherenow;
        prf_time tthen=tnow;
        wherenow=progr;
        tnow.update();
        rec_interval=(tnow-tthen);
        double current_speed=(wherenow-wherethen)/rec_interval;
        treq=(thegoal-wherenow)/current_speed;
        tavail=themaxtime-(tnow-time_strt);
        projected_pos=wherenow+((unsigned long) (current_speed*tavail));
        suspend_pos=wherenow+
                    ((unsigned long)(rec_interval*current_speed*
                                     ((int)(tavail/rec_interval))));

        if (tavail<rec_interval) advice=suspend;

        time_and_speed.push_back(std::make_pair(wherenow,current_speed));

        if (expecting_finish && tavail<treq) {
            report_asap=true;
            expecting_finish=false;
        } else if ((!expecting_finish)&& treq<tavail) {
            report_asap=true;
            expecting_finish=true;
        }
    }

    void Timer::forecast()
    {
        if (report_asap || wherenow >= milestones[ipt]) {
            if (thegoal<=wherenow) {
                prf::cout<<"Goal of "<<thegoal<<" has been reached.\n";
            } else {
                prf::cout<<posv<<" = "<<wherenow<<" after "
                <<Seconds(tnow-time_strt);

                if (limited) prf::cout<<". Remaining available time : "<<Seconds(tavail);

                if (treq<tavail) {
                    prf::cout<<". Estimated remaining time :"
                    <<Seconds(treq)<<"\n\n";
                } else {
                    prf::cout<<". Time needed to finish : "<<Seconds(treq)
                    <<".\nProjected final "<<posv<<": "<<projected_pos
                    <<". Estimated "<<posv<<" for suspend suggestion: "
                    <<suspend_pos<<"\n\n";
                }
            }

            if (wherenow>=(milestones[ipt]-mszone)&&ipt<npt) ++ipt;

            report_asap=false;

            prf::cout.flush();
        }
    }

    Timer::recommendation_type Timer::recommendation()
    {
        return advice;
    }

    void Timer::flush_timing_data(const char *speedfile)
    {
        prf::Output op(speedfile,"a");

        while (!time_and_speed.empty()) {
            op<<time_and_speed.front().first<<"\t"<<time_and_speed.front().second<<"\n";
            time_and_speed.pop_front();
        }

        op.close();
    }
}
