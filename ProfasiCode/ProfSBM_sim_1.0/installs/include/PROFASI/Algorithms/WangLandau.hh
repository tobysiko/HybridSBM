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

#ifndef WangLandau_HH
#define WangLandau_HH
#include "GMC.hh"
#include "../Observables/Observable.hh"

using std::vector;

namespace prf
{
    //! Wang Landau iterations
    /**
    This is an MC class derivative that redefines how the "Step" function
    works to implement the Wang Landau iterations. The acceptance probability
    for an update, which changes the energy from E1 to E2, is g(E1)/g(E2),
    where g(E) is an estimate of the density of states. The class keeps a
    histogram for this density of states and uses that in the Monte Carlo
    steps. A special feature of the Wang-Landau algorithm is that this
    estimate of the density of states is updated "live" during the simulation.
    The algorithm can be roughly described like this: \n\n
    <ol>
    <li> System is initialized randomly with an energy somewhere within the
         range</li>
    <li> A histogram of energy, H, is initialized in the range emin,emax </li>
    <li> Another histogram to store g(E) is initialized in the same range, but
        with every bin set to 1. </li>
    <li> A conformational update is performed</li>
    <li> Energies before and after the update are estimated</li>
    <li> Wang Landau acceptance criterion is applied, and the update is either
        accepted or rejected</li>
    <li> At the end of every update, g(Ei) for the bin corresponding to the energy
    of the system is multiplied by a certain factor, f. </li>
    <li> The appropriate bin of that histogram H is icremented by 1 </li>
    <li> The last 5 steps are repeated until a certain convergence criterion is met,
        defined as a certain degree of flatness of the occupation histogram. </li>
    <li> The factor f is reduced, and the occupation histogram (but not the
        estimate of g) is reset to zero</li>
    <li> Simulation continues from step 4, unless the factor f is already below
         a preset final value.</li>
    </ol>
    If the energy found is below the minimum in the range, the range is extended
    to accommodate the new value. Lower energy values are always interesting.
    If the update produces a new energy higher than the maximum, some heuristics
    is needed to bring it to within the range. This implementation can perform
    a series of high temperature Monte Carlo steps until the system re-enters
    the relevant range. This feature must be invoked with the switch --use_htmc.
    */

    class WangLandau : public GMC
    {
    public:
        WangLandau();
        ~WangLandau();
        //! Parse commands specific to the Wang Landau algorithm
        int parseCommand(InstructionString s);

        inline void energy_range(double e0, double e1) {emin=e0;emax=e1;}
        inline void flatness_cutoff(double v) { converr=v; }
        inline double flatness_cutoff() const { return converr; }
        inline double e_upper_limit() const { return emax; }
        inline double e_lower_limit() const { return emin; }
        inline void n_bins(int i) {nbins=i;}
        inline int n_bins() const { return nbins; }
        inline void init_g_factor() {gfact=gfact0;lngf=log(gfact);}
        inline void init_g_factor(double gf) { gfact=gfact0=gf; }
        inline double g(int ib) const {return gvals[ib];}
        inline double H(int ib) const {return his[ib];}
        inline double g_factor() const {return gfact;}
        inline void use_htmc(bool htm=true) {using_htmc=htm;}

        int update_g_factor();
        int check_bins();
        void init();
        //! Modify the MC step so that it perform Wang Landau updates
        int Step();
        int bin_of(double val);
        void init_g();
        void init_his();
        void reset_his();
        bool converged();
        int bring_to_range();
        int extend_till(double lowe);
        void save_state(std::string filename);
        void read_g(std::string filename);
        void print_setup();
        std::string ConfSignature();
    protected:
        double emin, emax, ebin, gfact, gfact0, lngf,converr;
        std::vector<double> gvals,his;
        int nbins,prevbins;
        bool using_htmc,usr_energy_range,usr_n_bins,gvalsread;
    };
}

#endif
