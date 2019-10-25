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

#ifndef HGMC_HH
#define HGMC_HH
#include "MC.hh"
#include <vector>

namespace prf {
    class HGMC : public MC
    {
    public:
        HGMC();
        virtual ~HGMC();
        Update *perform_update();
        virtual unsigned SwitchLambda();
        void set_n_temps(size_t i);
        void make_default_lambdas(size_t i);
        void set_lambda_range(double t0, double t1);
        //! Set temperatures from a given vector
        int SetLambdas(std::vector<double> tmpv);
        //! Specify temperatures in a file
        int SetLambdas(std::string tempfile);
        //! i'th beta value
        //inline double inverse_temperature(unsigned i) const {return beta[i%ntmp];}
        //! i'th temperature
        inline double lambda(unsigned i) const {return lambdas[i%ntmp];}

        //! Number of temperatures
        inline size_t number_of_lambdas() const { return ntmp; }
        //! Number of visits to the i'th temperature
        inline int NVisits(unsigned i) const {return dndT[i];}

        //! Index, or serial number of the current temperature
        inline unsigned CurLambdaIndex() const {return itmp;}

        //! Set current temperature
        inline void CurLambdaIndex(unsigned gind) {SetLambda(gind);}
        //! Set current temperature to the i'th
        void SetLambda(unsigned i);
         //! Write information about visits to different environments to file
        void writeLambdaStat(std::string flnm);
        void writeLambdas(std::string flnm);
        int prepare_continue(std::string flnm);
        //! Execute an instruction
        /**
          * The HGMC class handles some instructions of its own, as well as
          * forwarding some instructions to the UpdatesHandler. Instructions
          * for the HGMC class are:
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
        std::vector<double> lambdas;
        std::vector<size_t> dndT;
        unsigned ntmp,itmp,lambda_index;
        double tmin,tmax;
        bool tfileinuse,explicitntmps;
    };
}
#endif // HGMC_HH
