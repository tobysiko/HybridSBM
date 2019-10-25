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

#ifndef Bias_HH
#define Bias_HH
#include "Energy.hh"
#include <utility>

namespace prf
{
    //! Bias or E_{loc}
    /**
     * Bias is just a name for the energy term we have referred to as
     * $E_{loc}$.  It is a local electrostatic interaction between
     * adjacent peptide-bond-units.  In this implementation, the
     * interaction between the partial charges of the NH and CO
     * dipoles of two neighbouring peptide units are taken into
     * account. The term is not used for the two terminal residues,
     * unless they are protected by capping groups, called "EndGroup"s
     * in the program. A side-chain torsion potential and an
     * additional repulsion between neighbouring backbone H-H and O-O
     * pairs is also included. This is the simplest energy class in
     * PROFASI, and serves as a convenient template to write new
     * energy classes.
     */

    class Bias:public Energy
    {
    public:
        Bias();
        ~Bias();
        void init();
        double evaluate();
        double gradientXYZ(std::valarray<double> &gx);
        double deltaE(Update *);
        void Accept(Update *);
        void rangeEstimate(double &x1,double &x2);
        /// get & set
        inline void setHHRepulsionStrength(double k){ kappaH=k; }

        inline void setOORepulsionStrength(double k) { kappaO=k; }

        inline void setCoulombStrength(double k) { kcou=k; }

        inline double getHHRepulsionStrength() {return kappaH;}

        inline double getOORepulsionStrength() {return kappaO;}

        inline double getCoulombStrength() {return kcou;}

    private:
        double biasTerm(int chain, int residue);
        double repulsion(int chain, int residue);
        double repulsion_with_grd(int chain, int residue,
                                  std::valarray<double> &gx);
        double coulomb(int chain, int residue);
        double coulomb_with_grd(int chain, int residue,
                                std::valarray<double> &gx);
        double hrf(int i1,int i2,int i3,int i4, double strk);
        double dhrf(int i1,int i2,int i3,int i4,double strk,
                    std::valarray<double> &gx);
        double kappaH, kappaO, kcou, qC, qO, qN, qH;
        int nchanges;
        std::vector < std::vector < std::vector<double> > > q;
        std::vector < std::vector < std::vector<int> > > qPos;
        std::vector < std::vector<double> > Vebias, Mebias;
        std::vector< std::pair<int, int> > changed;
    };
}

#endif
