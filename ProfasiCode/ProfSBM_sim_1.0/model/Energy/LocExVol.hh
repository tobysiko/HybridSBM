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

#ifndef LocExVol_HH
#define LocExVol_HH
#include "Energy.hh"
#include "ExVolBase.hh"

namespace prf
{
//! Third-neighbour excluded volume contribution
    /**
     * \ingroup profasi_energies
     * Contribution to excluded volume from "all" atom pairs separated by three
     * covalent bonds. "All" excludes pairs with fixed distance, like in phenyl
     * rings. All such pairs are listed during initialization, and calculation
     * proceeds by a straight forward summation over such pairs. During an update
     * of a single angle, at most 9 such pairs are affected. What pairs would be
     * affected by different kinds of updates is also pre-calculated and stored.
     * Very little needs to be done to get deltaE after an update.
     *
     * In an equilibrated system this term contributes the lion's share of
     * excluded volume energy. Because the atom pairs involved in this term can not
     * fly apart. This term has the biggest influence on the Ramachandran angle
     * distributions and also the side chain chi angle distributions.
     *
     * By definition this term always involves atoms from the same peptide. So,
     * it does not need to use the periodic distance measure. By design, the size
     * of the periodic box is sufficiently large so that the different parts of
     * a chain can not reach out and interact with their copies in a different
     * periodic box. And atoms of one chain are always kept together, so that
     * separation between them can not exceed the box length.
     *
     */

    class LocExVol:public ExVolBase, public Energy
    {
    public:
        LocExVol();
        ~LocExVol();
        void init();
        void rangeEstimate(double &x1, double &x2);
        double evaluate();
        double gradientXYZ(std::valarray<double> &ans);
        double deltaE(Update *);
        void Accept(Update *);
        bool is_loc_pair(int i, int j);
        inline int NPairs() const {
            return npair;
        }

        void PrintPairs() const;
    private:
        std::valarray < double >sig2, asa, bsa;
        std::valarray < double >Melpsa, Velpsa;
        std::valarray < int >lci1, lci2, lcid;
        std::valarray < int >doftorng;
        int npair, maxdof, totdof, upairbegin, upairend;
        double Vexv(int i);
        double dVexv(int i, std::valarray<double> &gx);
    };
}

#endif
