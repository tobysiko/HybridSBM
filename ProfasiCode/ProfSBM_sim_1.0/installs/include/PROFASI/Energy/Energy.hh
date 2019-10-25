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

#ifndef Energy_HH
#define Energy_HH
#include "../Elements/Population.hh"
#include "../Updates/Update.hh"
#include "../Observables/Observable.hh"

/**
* \defgroup profasi_energies Energy Classes
*/

namespace prf
{
    //! Energy base class
    /**
     * The Energy base class declares an interface for different energy terms
     * which can be used, for instance, to decide whether a certain update
     * is to be accepted. All energy terms in ProFASi are implemented as
     * derived classes of this class.
     * \ingroup profasi_energies
     */

    class Energy : public Observable
    {
    public:
        virtual ~ Energy();
        //! Connect the energy term with one population
        void Connect(Population * pl);
        //! Set some adjustible parameters
        /**
          * This function provides a unified mechanism to adjust parameters
          * in an energy term. What parameters it might make sense to adjust
          * depends on the energy term. If the energy term is a pseudo-energy
          * implementing RMSD constraints relative to a structure abc.pdb,
          * the name of the structure file should be an easily adjustible
          * parameter rather than a hard coded string. The H--O distance
          * in a hydrogen bond, on the other hand, is less useful as an
          * adjustible parameter. It is up to the individual derived classes
          * to do what they like with the requested parameter changes. The
          * base class ignores the parameter requests.
          */
        virtual void set_pars(std::string pars);
        virtual void init();
        // In case you wish to initialize, although that was done before
        void re_init();
        //! Overrides the virtual member from the Observable class.
        void refresh();

        //! return result from last calculation
        inline double value() {
            return vval;
        }

        //! return result from the last energy change calculation
        inline double deltaE() {
            return delv;
        }

        //! do a new ab initio calculation
        virtual double evaluate(); //!<do a new ab initio calculation
        //! do a new delta calculation optimized for the given kind of update
        virtual double deltaE(Update *);
        double delta(Update *);
        //! Calculate energy change for an update, but with a stop condition
        /**
          * This function contains one of the interesting speed-up tricks of
          * ProFASi. In Markov chain Monte Carlo simulations, there is an
          * acceptance criterion, which depends, on energy changes, temperature
          * a random number etc. Often, it is possible to calculate everything
          * other than the energy change before proceeding to the energy
          * change calculations. This means, before starting out to find how
          * much the energy changed, we can figure out the maximum acceptable
          * total energy change. We then isolate the most expensive contribution
          * and calculate the regular delta E for all other terms, and subtract
          * from the maximum total acceptable change, and get a value for the
          * maximum allowed value for the last and most expensive term. We pass
          * this value to the delta E calculation of the last term, using
          * this function. Implementation of that energy function is allowed
          * to abort its calculation, if at some point during its calculation,
          * it can determine that there is no way to stay below the given limit.
          * This means, this function is allowed to return any (wrong) value
          * greater than maxde, if the correct value is also greater than maxde
          * with certainty. It is understood that the update will be rejected
          * if a value greater than maxde is returned.
          *
          * As an example, consider the situation when the most
          * expensive term is positive definite. We start by calculating the
          * contributions of the moved atoms before the update. Then while
          * calculating the contributions of the same atoms after the update, if
          * after 5 of the 2000 terms, the partial sum exceeds a limit, there is
          * no point in calculating the remaining contributions, which can not
          * bring the total change below maxde. This trick makes a noticeable
          * impact on the execution speed.
          */
        virtual double deltaEwithlimit(Update *updt, double maxde);
        //!< Like deltaE, but stop if change exceeds a given limit

        //! Accept a proposed update
        /**
         * PROFASI uses many local contribution matrices and backup
        * variables to quickly evaluate the change caused to any energy
        * term by a given update. Those backup variables and the
        * contribution matrices need to be kept in sync with the current
        * system configuration. So, each energy class has an accept and
        * a reject method, which are called when an update is accepted
        * or rejected. What is done inside these functions depends on
        * the energy term and optimization trick used therein.
         */
        virtual void Accept(Update *);
        //! Reject a proposed update
        virtual void Revert(Update *);
        inline int NC() const {
            return p->NumberOfChains();
        }

        virtual void rangeEstimate(double &x1,double &x2);
        bool is_in_use() const {return inuse;}
        inline void is_in_use(bool flg) { inuse=flg; }
        int init_obs();
    protected:
        Energy();
        double vval, delv;
        bool inuse,initialized;
    };
}

#endif
