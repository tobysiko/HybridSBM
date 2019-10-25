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

#ifndef TorsionTerm_HH
#define TorsionTerm_HH
#include "Energy.hh"
#include <utility>
#include <stack>

namespace prf
{
    //! Torsion angle potential
    /**
     * Distributions of various torsional degrees of freedom in PROFASI
    * is approximately correct without the need of a term in the energy
    * that explicitly depends on the torsion angles. The purely repulsive
    * excluded volume term and the local backbone electrostatics lead to
    * reasonable distributions of the torsion angles. But there are a few
    * that misbehave. The angular distributions of Asparagine, the backbone
    * angle distribution of Glycine are important examples. It indicates
    * that there is some relevant physics for a few amino acids, which is
    * missed by a model based purely on repulsions. This term is to address
    * only those cases. It is to be interpreted as a proxy for some physics
    * we do not understand perfectly.
    * \ingroup profasi_energies
    */

    class TorsionTerm: public Energy
    {
    public:
        TorsionTerm();
        ~TorsionTerm();
        void init();
        double evaluate();
        double evaluate(Update *);
        double gradientDOF(std::valarray<double> &ans,
                        std::vector<int> &indxs);
        double deltaE(Update *);
        void rangeEstimate(double &x1,double &x2);
    private:
        std::valarray<int> location,torstype,dofid;
        void regtors(int i, int &itors,int ityp);
        double torsion(int dofsno);
        double dtorsion(int dofsno);
        double torFcn(double phi, int type);
        double dTorFcn(double phi, int type);
        double kgly;
        int ntors;
    };
}

#endif
