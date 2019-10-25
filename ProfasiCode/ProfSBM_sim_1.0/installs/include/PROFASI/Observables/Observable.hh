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

#ifndef Observable_HH
#define Observable_HH
#include "../Aux/Named.hh"
#include <vector>
#include "../Aux/profasi_io.hh"
#include "../Aux/AdaptiveHis.hh"
#include "../Aux/statistics.hh"
#include "../Updates/Update.hh"

/**
* \defgroup profasi_observables Measurements
* Measurement of different properties of the protein chain or more general
* properties of the system are performed through "Observable" classes. An
* Observable typically calculates one specific property of the system, and
* returns a double precision value when asked.
*
* All Observables in PROFASI (version 1.1 and later) can be initialized
* and customized through a string based interface. They can receive "commands"
* for specific settings as strings, and store them. At the time of
* initialisation, they make use of the commands they have received to
* customize their properties. Of course, the commands have to make sense for
* an Observable to implement them. If you send and arbitrary command, such as
* "set cutoff -56.7" to an observable that calculates the radius of gyration,
* it will simply be ignored. But something like "histogram -10.0 10.0 100"
* might potentially make sense for most measurements.
*/

namespace prf
{
    //! An observable is in principle anything that is named and has a value
    /**
     * Common interface for a variety of "observables". In PROFASI, an
     * Observable is simply anything that has a name and a value. All kinds
     * of energies, RMSD, secondary structure content ... may fall in this
     * category. While the method by which each individual observable is
     * calculated may vary, as well as its role in the program as a whole,
     * for an outside analysis, such details often don't matter. So, it is
     * often desired to make an array of observables and do one common
     * operation on them. This would not be possible if, for instance, Energy
     * and ProteinRMSD were incompatible types. But since both inherit from
     * this observable class, they are compatible classes as far as the
     * properties defined here are concerned. If you wish to implement a new
     * observable and use it like other Observables of PROFASI, just inherit
     * from this class, and make sure you implement a refresh() function.
     * The refresh function should calculate the value and assign it to a
     * member variable called obsval. Sometimes, the value would be calculated
     * in some other context. But nevertheless, it is necessary to implement
     * this refresh function, even if all it does is to assign the already
     * calculated value to obsval.
     * \ingroup profasi_observables
     */

    class Observable : public Named
    {
    public:
        Observable();
        virtual ~Observable();
        //! Give the object an instruction to process during initialization
        void set(std::string cmd);
        //! All observables must implement one initialize routine
        /**
        * Even if it seems that one particular observable might need additional
        * arguments during initialization, it is advantageous to have a uniform
        * syntax for all of them. So, when additional arguments are needed,
        * one should provide them in a separate function called before
        * initialization, and then call init_obs without arguments. The name
        * init_obs instead of a more natural "init" or "initialize" is because
        * an Observable often inherits from other classes which represent its
        * character more fundamentally. So, the names such as "init" are
        * kept free for such base classes.
        */
        virtual int init_obs();

        //! Necessary before an observable value is used.
        /**
         * This is done because complex observables like RMSD are not evaluated
         * at every step. A call to refresh() would make sure that the
         * Observable has its most current value. The optional argument tindx
         * was introduced in version 1.1 when management of histograms was
         * relocated from ObsHandler class to the Observable class. The argument
         * tells the Observable object about a "temperature index" which it can
         * use to put the current data in the appropriate histogram block.
         */
        virtual void refresh();
        //! Actual calculation of the value of the observable
        virtual double evaluate();
        virtual double gradientDOF(std::valarray<double> &ans,
                                std::vector<int> &indxs);
        virtual double gradientXYZ(std::valarray<double> &ans);
        //! Quick estimate of the change in an Observable due to an update
        /**
          Estimate how much the value will change due to the a given update.
          Quite often it is possible to estimate the change in the value of
          an Observable due to such an MC update at a much smaller computational
          cost than evaluating the Observable from scratch.

          This function assumes that the observable was evaluated at one
          state of the population, and the update passed as argument was
          performed on that state of the population. It is also ok if the
          state of the population is only changed with MC updates and the
          accept or reject function is called for each accepted or rejected
          update.

          Except for energy classes in ProFASi, the other Observables are not
          evaluated or kept up to date during an MC sweep. For such variables
          this delta function only returns something useful immediately after
          an evaluate call.
          */
        virtual double delta(Update *u);
        virtual void accept(Update *u);
        virtual void reject(Update *u);
        //! Estimate a range in which values of this observable are expected
        /**
         * The default is between 0 and 1. So, for observables with values
         * always between 0 and 1, you need not over-write this virtual
         * function. Sometimes the observable will have a different fixed
         * range, determined by its definition. Sometimes the range can not
         * be determined perfectly. In such a case, let this function just
         * return something reasonable.
         */
        virtual void rangeEstimate(double &x1,double &x2);
        //! Retrieve the value of the observable
        inline double operator()() {return obsval;}

        //! Retrieve the value of the observable
        inline double Value() {return obsval;}

        inline prf_utils::AdaptiveHis * histogram() { return myhis; }

        inline void output_prefix(std::string prx) {oprefx=prx;}

        inline std::string output_prefix() const {return oprefx;}

        virtual void write_rtkey(Output &op);
        virtual void write_snapshot(Output &op);
        virtual void his_range(double xmn, double xmx);
        virtual void his_bin_size(double sz);
        inline double his_bin_size() const {return xbin0;}
        virtual void his_setup();
        virtual void his_fill(int i);
        virtual void his_reset();
        //! Stop collecting statistical data like averages and histograms
        inline void disable_stats() { gathstat=false; }
        //! Start collecting statistical data like averages and histograms
        inline void enable_stats() { gathstat=true; }
        inline void make_his(bool sw) {nohis=(not sw);}
        inline void set_his_resume(bool sw=true) { hisresume=sw; }
        virtual void his_save();
//        virtual void his_adjust();
        virtual void avg_fill(int i);
        virtual void avg_reset();
        virtual void avg_write(Output &op);
        void set_n_temp(int i);

        inline void setPopulation(Population *popl) {p=popl;}
        inline int grad_eval_type() const { return grdtyp; }
        virtual void set_logger_threshold(int thr);

    protected:
        double obsval, obsdel;
        std::vector<std::string> usrcmd;
        bool requires_population, fixed_his, gathstat, nohis;
        bool userhisrange,usernbins,userbinsz,hisresume;
        Population *p;
        int ntmp,nbins0,grdtyp;
        double xmin0,xmax0,xbin0;
        std::vector<prf_utils::StatBox> avgs;
        prf_utils::AdaptiveHis *myhis;
        std::string oprefx;
        int auto_nbins(int nexpected);
        int pop_check();
        int log_thres;
    };
}

#endif
