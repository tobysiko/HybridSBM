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

#ifndef GMC_HH
#define GMC_HH
#include "MC.hh"
#include <vector>

namespace prf {
    class GMC : public MC
    {
    public:
        GMC();
        virtual ~GMC();
        Update *perform_update();
        virtual unsigned SwitchTemp();
        void set_n_temps(size_t i);
        void make_default_temps(size_t i);
        void set_temperature_range(double t0, double t1);
        //! Set temperatures from a given vector
        int SetTemps(std::vector<double> tmpv);
        //! Specify temperatures in a file
        int SetTemps(std::string tempfile);
        //! i'th beta value
        inline double inverse_temperature(unsigned i) const {return beta[i%ntmp];}
        //! i'th temperature
        inline double temperature(unsigned i) const {return 1.0/beta[i%ntmp];}

        //! Number of temperatures
        inline size_t number_of_temperatures() const { return ntmp; }
        //! Number of visits to the i'th temperature
        inline int NVisits(unsigned i) const {return dndT[i];}

        //! Index, or serial number of the current temperature
        inline unsigned CurTempIndex() const {return itmp;}

        //! Set current temperature
        inline void CurTempIndex(unsigned gind) {SetTemp(gind);}
        //! Set current temperature to the i'th
        void SetTemp(unsigned i);
         //! Write information about visits to different environments to file
        void writeTempStat(std::string flnm);
        void writeTemperatures(std::string flnm);
        int prepare_continue(std::string flnm);
        //! Execute an instruction
        /**
          * The GMC class handles some instructions of its own, as well as
          * forwarding some instructions to the UpdatesHandler. Instructions
          * for the GMC class are:
          * \li \b --temperature_file or \b -tfile : If you want temperatures
          * to be read in from a certain file.
          * \li \b --num_temperatures or \b -ntmp :  Number of temperatures
          * Example: -ntmp 10
          * \li \b --min_temperature or \b -tmin : Minimum temperature
          * Examples: -tmin 0.42 or -tmin "280 Kelvin"
          * \li \b --max_temperature or \b -tmax : Maximum temperature
          */
        virtual int parseCommand(InstructionString s);
        virtual void print_setup();
        virtual std::string ConfSignature();
    protected:
        std::vector<double> beta;
        std::vector<size_t> dndT;
        unsigned ntmp,itmp;
        double tmin,tmax;
        bool tfileinuse,explicitntmps;
    };
}
#endif // GMC_HH
