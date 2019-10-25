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

#ifndef Hydrophobicity_HH
#define Hydrophobicity_HH
#include "Energy.hh"

namespace prf
{
//! Effective hydrophobic attraction between non-polar side chains
    /**
     * This class represents the hydrophobicity interaction in the model. In its
     * current version, hydrophobicity is a simple pair-wise additive potential.
     * A degree of contact between two hydrophobic side chains is calculated
     * based on the proximity of a pre-determined set of atoms from each
     * side chain. This degree of contact is multiplied by a strength assigned
     * to pairs of side chain, which are in rough correspondence with their size.
     * Like other terms, the implementation of the class takes care of the
     * optimization. Contribution from each pair of hydrophobic side-chains is
     * saved in a contribution matrix. For each update, the affected
     * contributions are determined, and reevaluated, and compared with the
     * value stored in the matrix. Contributions from hydrophobic pairs that
     * are nearest and second-nearest neighbours in sequence are supressed by
     * factors 0.0 and 0.5 respectively.
     * \ingroup profasi_energies
     */

    class Hydrophobicity:public Energy
    {
    public:
        Hydrophobicity();
        ~Hydrophobicity();
        void init();
        double evaluate();
        double gradientXYZ(std::valarray<double> &gx);
        double deltaE(Update *);
        void Accept(Update *);
        void DisplayMatrix();
        void rangeEstimate(double &x1, double &x2);
        inline void SetNNStrength(double n1, double n2) {
            strnn = n1;
            strnnn = n2;
        }

        inline void ScaleStrength(double x) {
            sclfct = fabs(x);
        }

        inline double contribution(int i, int j) const {
            return sclfct*Mehp[NHAA * i + j];
        }

        inline int nhaa() const {
            return NHAA;
        }
        inline void setfocalRes(int res){focalRes=res;}
        inline void setfocalRes2(int res){focalRes2=res;}
        inline void setfocalResScale(double factor){focalResScale=factor;}
        
        double interhp();
        //! i and j refer to the list of hydrophobic residues in the whole system.
        double contact_frac(int i, int j,double(*distf)(int,int));
        //! iaa and jaa refer to the ligand indices relative to the whole system
        double hp_contact_frac(int iaa, int jaa,double(*distf)(int,int));
        double InterChain(int ich, int jch);
    private:
    	int focalRes,focalRes2;
    	double focalResScale;
        int htyp(OneLetterCode cod);
        int groupOf(OneLetterCode cod);
        double hp_pair(int i, int j,double(*distf)(int,int));
        double dhp_pair(int i, int j, bool perhint,
                        std::valarray<double> &gx);
        std::valarray < double > Mehp, Vehp;
        std::valarray<int> changed;
        int NHAA, N, nchanges;
        double a2, b2, d2,slpe;
        double strnn, strnnn, sclfct,max_at_dist2,max_cmp_dst2;
        std::valarray<double> strength;
        std::valarray<int> atom;
        std::valarray<int> natom, ih, chainof, soften,hlig;
        double r2min[12];
        int r2minpartner[12];
        std::valarray <OneLetterCode> seq;
        static const int nHPgrp;
        static const double hpstr[];
        static std::string hatms[][6];
    };
}

#endif
