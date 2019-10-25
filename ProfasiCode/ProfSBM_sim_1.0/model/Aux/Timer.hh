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

#ifndef Timer_HH
#define Timer_HH
#include <cstdio>
#include "prf_time.hh"
#include <cstdlib>
#include <string>
#include <list>
#include <vector>

namespace prf_utils
{
    //! ProFASi's own timer for managing time limited runs
    /**
      * Meaningful simulations of interesting protein systems take a long
      * time. Even with the simplified force field and optimized code of
      * ProFASi, a run might require weeks of computing time. But on large
      * computing facilities, there are limits to how long a single job
      * can run, let's say 12 hours.
      *
      * ProFASi's simulation programs, such as ParTempRun or SimTempRun
      * make use of this timer class to cleanly suspend an on going run
      * when the allocated time is about to finish. The runs are started
      * with a certain amount of "available time", passed to the programs
      * through the settings file or the command line arguments. The timer
      * measures how fast the program is going about its task. It estimates
      * how far it can go in the available time. If it estimates that it
      * will not come back to the next "check point" (configuration save)
      * before the allocated time is up and the run is killed by the
      * system, it recommends that the simulation program "suspends".
      *\ingroup utilities
      */

    class Timer
    {
    public:
        //! Timer recommendation type
        typedef enum {go_on, suspend} recommendation_type ;
        //! Default constructor
        Timer();
        ~Timer();
        //! Start the stop-watch
        void Start(unsigned long startpos=0);
        //! Stop the stop-watch
        void Stop(unsigned long endpos=0);
        //! Reset the stop-watch
        void Reset();
        //! Interval since the start or the last reset in seconds
        double Measurement();
        //! Print interval between start and stop in hours, minutes, seconds etc.
        void PrintInterval();
        //! A name for what is changing with time
        /**
          * If we are going somewhere, we might estimate how far we have come and how
          * much of the journey remains in kilometres. For an MC simulation, the progress
          * is measured in MC cycles. This function puts a label on the space axis.
          *
          */
        inline void space_axis_name(std::string nam) {posv=nam;}

        inline std::string space_axis_name() const {return posv;}

        //! Set a new target to be reached
        inline void set_new_goal(unsigned long newgoal) {thegoal=newgoal;}

        //! Tell the timer how much progress has been made
        inline void set_current_progress(unsigned long progr) {wherenow=progr;}

        //! Tell the timer that we have infinite time for the task
        /**
          * If there is unlimited available time for the task, the timer does not
          * calculate whether or not the job should be suspended. It merely
          * reports how much longer the job might take, when asked.
          */
        inline void have_endless_time() { limited=false; }

        //! Tell the timer how much time is available for the task
        inline void allocate_time(double dur) { themaxtime=dur; limited=true; }

        //! Write data about how much time was spent in different stages of the task
        void flush_timing_data(const char *timefile);

        //! Set how often the timer should talk about the progress made, time left etc.
        inline void n_progress_reports(int ntpt) {npt=ntpt;}

        //! Record progress at intervals of progress
        /**
          * At regular intervals during the execution of the task, this function
          * should be called with a value indicating how much of the task has
          * been accomplished. The timer has no way of knowing how many MC cycles
          * have been finished or how many km have been travelled. But the
          * application can every once in a while tell the timer what has happened.
          * The timer will associate a time with that, and use the information
          * to calculate what might happen in the remaining available time.
          */
        void record(unsigned long progr);
        //! Calculate how much more time is needed for the job etc.
        /**
          * This function calculates, how much more time is needed for the job,
          * and if the job can not be completed in the remaining available time,
          * how much more work can be expected to be done in the remaining time.
          */
        void forecast();
        //! Ask the timer, if the application should continue or stop
        recommendation_type recommendation();
        //! Has the goal been reached ?
        inline bool goal_reached() {return wherenow>=thegoal;}

        //! Tell the timer that we are not going any further
        /**
          * Sometimes, although we expect it would take 1 billion MC cycles
          * to do something, the simulation program might decide that
          * enough is enough, and we have found a converged answer. It would
          * then abandon the mission of reaching the 1 billion cycles. If
          * it tells the timer that the mission is abandonned, the timer
          * prints a more fitting closing message.
          */
        inline void abandon() { mission_abandonned=true; }

    private:
        prf_time time_strt, time_tmp, time_end, tnow;
        unsigned long thegoal, thestart, wherenow, projected_pos,
        suspend_pos;
        double themaxtime,treq,tavail, rec_interval, mszone;
        int npt,ipt;
        bool limited, report_asap, expecting_finish,mission_abandonned;
        recommendation_type advice;
        std::string posv;
        std::list<std::pair<unsigned long, double> > time_and_speed;
        std::vector<unsigned long> milestones;
    };
}


#endif
