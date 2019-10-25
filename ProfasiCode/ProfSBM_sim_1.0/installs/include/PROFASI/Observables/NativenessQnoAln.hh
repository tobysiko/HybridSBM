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

#ifndef NativenessQnoAln_HH
#define NativenessQnoAln_HH
#include "Observable.hh"
#include "../Elements/Population.hh"
#include "../Aux/PDBReader.hh"

namespace prf
{
    //! A measure of nativeness based on local environment of residues
    /**
    * This class implements the nativeness parameter Q introduced in
    * S. Takada, Z.A. Luthey-Schulten, P.G. Wolynes, J. Chem. Phys.
    * 110 (1999) 11616. The definition of the local version of the
    * measure is taken from G. Chikenji, Y. Fujitsuka and S. Takada, Chem.
    * Phys. 307 (2004) 157.
    *
    * This measure is only concerned with Calpha distances. For each
    * residue, a list is made of all the other residues which are within
    * 6 Angstroms CA distance from it in the native state. To assess the
    * nativeness of a given new state, each residue's new distance from
    * those in its native neighbourhood list is evaluated. For each pair,
    * the quantity exp(-(r_{ij}-r0_{ij})^2 / d^2) is calculated. These
    * are then averaged over all neighbours of a residue to define its
    * local "nativeness": how similar its environment in the current state
    * is, compared to the native state. The average over all unique pairs
    * over the entire protein is the value returned by the observable.
    * \ingroup profasi_observables
    *
    */

    class NativenessQnoAln : public Observable
    {
    public:
        NativenessQnoAln();
        ~NativenessQnoAln();
        int init_obs();

        inline void nativeStructure(std::string pdbselection)
        {infile=pdbselection;}

        int init();
        double evaluate();
        inline int n_loc_vals() const {return (int) res.size();}

        void loc_val(size_t i, int & resi, double &qi);
        void avg_fill(int i);
        void avg_reset();
        void avg_write(Output &op);
        // Overwrite base class function to not let the range change
        void his_range(double xmn, double xmx);
    private:
        std::string infile,altcrd;
        std::vector<std::vector<int> > neighbour;
        std::vector<std::vector<double> > refdist;
        std::vector<int> res, n_neighbours, cas, nevals;
        std::vector<double> Qloc;
        Matrix<double> Qlochist;
        std::valarray<unsigned long> dndt;
        int ntotpairs,nres;
        double d0,dcut,scale2;
    };
}

#endif
