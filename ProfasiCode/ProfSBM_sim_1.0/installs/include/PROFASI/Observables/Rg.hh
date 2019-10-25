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

#ifndef RG_HH
#define RG_HH
#include "Observable.hh"
#include "../Elements/Population.hh"

namespace prf
{
    //! Radius of gyration of a section of a protein chain
    /**
    * This observable returns the radius of gyration of all non-hydrogen
    * atoms of a given chain or of the whole population. In the case of
    * rg over many chains, it is assumed that EnforceBC has been called,
    * so that all coordinates are within one unit periodic box. Only
    * non-hydrogen atoms are used in the calculation.
    *
    * \ingroup profasi_observables
    */

    class Rg : public Observable
    {

    public:
        Rg();
        ~Rg();
        int init_obs();

        //! Limit to chain i. If i is -1, all chains.
        inline void of_chain(int i) {ich=i;}

        double evaluate();
        double gradientXYZ(std::valarray<double> &ans);
        void rangeEstimate(double &x1, double &x2);
    private:
        int ich;
    };
}

#endif

