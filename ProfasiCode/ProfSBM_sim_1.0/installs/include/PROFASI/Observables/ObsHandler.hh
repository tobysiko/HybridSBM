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

#ifndef ObsHandler_HH
#define ObsHandler_HH
#include "Observable.hh"
#include "../Energy/Extras/ObsEnergy.hh"
#include "../Aux/His2D.hh"
#include "../Aux/HandlerBase.hh"
#include <vector>
#include <map>
#include <string>
#include <algorithm>

using namespace prf;

using namespace prf_utils;

namespace prf
{
    //! A handler class for Observables
    /**
     * The ObsHandler, or observable handler takes care of routine tasks in
     * connection with taking measurements and performing three frequent
     * operations on them.
     */

    class ObsHandler : public HandlerBase
    {

    public:
        ObsHandler();
        ~ObsHandler();
        void setPrefix(std::string file_prefix);
        void setDir(std::string write_dir);
        void track(Observable *ob,std::string options);
        int set_block_props(int num_states,std::string state_designator);
        int initialize();
        int parseCommand(InstructionString s);
        void make_obs(InstructionString spec);
        void set_obs_props(InstructionString spec);
        int add_2d_his(InstructionString s);
        void printConfig();
        void sample(int current_time,int current_state);
        void writeRTKey();
        void writeRTSnapshot();
        void writeHistograms();
        void writeAverages();
        void disableStatistics();
        void enableStatistics();
        void resetHistograms();
        void resetAvgData();
        void restoreHistograms(std::string filename);
        void restoreAvgData(std::string filename);
        void saveHistograms();
        void saveAvgData();

        inline int num_obs() const {return nobs;}

        inline void setPopulation(Population *pl) { popl=pl; }

        inline void setSwitch(std::string sw, bool val=true)
        {swtch[sw]=val;usrspc[sw]=true;}

        inline void unsetSwitch(std::string sw) {setSwitch(sw,false);}

        int make2DHis(std::string obsname1, std::string obsname2);
        inline His2D *get2DHis(int i,int itmp) {return hs2d[i*nstates+itmp];}

        void setup2DHis();
        void save2DHis();
        void save2DHis(int i,int j);
        Observable * get_obs(std::string nm);
        Observable * get_obs(int i) { return (i>=0&&i<nobs)?obs[i]:NULL; }
        void re_init_obs();
        void fix_dependencies(ObsEnergy *depterm);

    private:
        void add_obs_and_opts(Observable *tmp, std::string nm, std::string pars);
        void use_obs(Observable *o, bool flg=true);
        std::string myprefx,avgfile,rtfile,dirc,stdesg;
        std::map<std::string,int> obsno;
        std::map<std::string,bool> swtch,usrspc;
        std::vector<Observable *> obs, unused, hs2dobs1,hs2dobs2,myobs;
        std::vector<std::string> hs2dreq1,hs2dreq2;
        int nstates,curst,curtime,nobs;
        std::vector<His2D *> hs2d;
        Population *popl;
        Output rt;
    };

}

#endif
