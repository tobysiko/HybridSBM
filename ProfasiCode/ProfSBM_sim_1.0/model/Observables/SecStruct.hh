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

#ifndef SecStruct_HH
#define SecStruct_HH
#include "Observable.hh"
#include "../Elements/Population.hh"
#include <vector>

namespace prf
{
    //! Secondary structure analyzer
    /**
    This is a simple but extremely versatile class to analyze secondary
    structure based on Ramachandran angles. Essentially, it has a list of
    residues whose behaviour it tracks, and during a refresh, it checks
    which of the residues are in a certain box (phimin to phimax, psimin
    to psimax) in the Ramachandran plot. The boundaries of the box can be
    set by the user. One can choose "default_limits helix" or
    "default_limits strand" to choose the PROFASI defaults to regard a
    residue as a helix or a strand. One can pass any set of residues to
    track in a file. The file should contain the amino acid indices of
    the relevant residues relative to the entire population. One can
    instruct the Observable to use the residues of one chain. If the user
    makes no choices about the residues, the entire population is chosen.
    But, in these last two cases, the amino acids at the ends of each
    chain are omitted.

    This single class can be used to construct Observables to measure helix
    content, beta strand content, native helix content, native strand
    content, and even something to monitor a totally new area of the
    Ramachandran plot for the frequency of visits by different residues.

    For every one of the uses mentioned above, the class also creates a
    "profile" curve, where residue-wise average properties are saved as
    a function of temperature index in a file with the suffix ".profile".
     */

    class RCBin : public Observable
    {
    public:
        RCBin();
        ~RCBin();
        inline void set_limits(double phmin,double phmax,
                               double psmin,double psmax) {
            phimin=phmin;phimax=phmax;psimin=psmin;psimax=psmax;
        }

        int set_res_list(std::string filename);
        double evaluate();
        int init_obs();
        void avg_fill(int i);
        void avg_reset();
        void avg_write(Output &op);
        inline int is_in_bin(int i) { return sttres[i]; }

        // Overwrite base class function to not let the range change
        void his_range(double xmn, double xmx);
        void his_nbins(int n);
        void rangeEstimate(double &x0, double &x1);
    private:
        double phimin,phimax,psimin,psimax;
        std::valarray<unsigned long> dndt;
        std::vector<int> refvec, sttres;
        Matrix<double> prof;
    };
}

#endif
