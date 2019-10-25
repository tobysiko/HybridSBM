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

#ifndef OligoOrient_HH
#define OligoOrient_HH
#include "Observable.hh"
#include "SecStruct.hh"
#include "../Elements/Population.hh"
#include "ContactFunctions.hh"

namespace prf
{
    //! Number of pairs of parallelly oriented peptides in oligomers
    /**
     * This Observable counts the number of pairs of peptides in the system
     * with at least 3 inter-chain backbone-backbone hydrogen bonds and xx%
     * (adjustible) beta strand content in both, which are oriented parallel
     * to each other. Parallel, here, means that the end-to-end vectors of the
     * peptides have an angle less than dd (adjustible) degrees.
     *
     * The Observable makes no sense for single peptide runs, and it may not
     * make sense to use it for every oligomerisation simulation. It was
     * originally made for the study of the 7 residue peptide Abeta(16-22),
     * but for longer systems, the alignment or otherwise of end-to-end vectors
     * might mean much less.
     *
     * Normally one is interested in the content of both parallel and
     * anti-parallel pairs. So, this class actually finds both. But as an
     * observable it only returns the number of parallel pairs in the function
     * Value(). But value2() returns the anti-parallel pair content. To create
     * an Observable object out of the number of anti-parallel pairs,
     * a dummy class OligoNminus is provided below, which depends on this
     * class to get its value.
     * \ingroup profasi_observables
     */

    class OligoOrient : public Observable
    {
    public:
        OligoOrient();
        ~OligoOrient();

        double evaluate();
        int init_obs();
        void his_fill(int curT);
        void avg_fill(int curT);
        void his_reset();
        void avg_reset();
        void his_save();
        void avg_write(Output &op);
        void write_rtkey(Output &op);
        void write_snapshot(Output &op);
        //! Angle cut-off below which to strands should be regarded as parallel
        /**
         * Should be liberal. The default is 30 degrees. Values are to be
         * passed in radians, not degrees.
         */
        inline void setAngleCut(double x) {copang=cos(x);}

        //! Minimum energy of backbone-backbone hydrogen bonds
        /**
         * Minimum energy of inter-chain backbone-backbone hydrogen bonds
         * for a pair to be considered bonded in a dimer, or as a part of
         * a larger oligomer. The default is 4.6, corresponding to 2-3
         * hydrogen bonds.
         */
        inline void minHB(double mbnd) {hbcut=mbnd;}

        //! Minimum beta strand content for each peptide
        /**
         * Peptides with less beta strand content than what is provided here
         * will not be considered in the calculations. If nothing is provided,
         * the default value is 0.5.
         */
        inline void minBeta(double mbt) {minbt=mbt;}

        void rangeEstimate(double & x1, double &x2);
        void his_setup();
    private:
        HBContactChains cf;
        std::vector<RCBin> chbeta;
        int np,nm;
        double copang,minbt,hbcut;
        std::vector<prf_utils::StatBox> avgb;
        prf_utils::His1D *hisb;
    };

}

#endif
