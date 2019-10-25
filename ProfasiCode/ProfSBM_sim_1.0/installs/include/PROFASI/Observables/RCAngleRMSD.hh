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

#ifndef RCAngleRMSD_HH
#define RCAngleRMSD_HH
#include "Observable.hh"
#include "../Elements/Population.hh"

namespace prf
{
    //! Root mean square deviation of Ramachandran angles
    /**
    * This is a similarity measure based on the backbone degrees of
    * freedom in the model, the Ramachandran phi and psi angles. The
    * return value of the observable is the root mean square difference
    * from a given reference set of angles. Residuewise differences
    * can also be retrieved. The angular differences are of course,
    * modulo pi.
    *
    * The reference phi and psi angles for each residue can be given
    * in a file, or calculated from a PDB file.
    * \ingroup profasi_observables
    */

    class RCAngleRMSD : public Observable
    {
    public:
        RCAngleRMSD();
        ~RCAngleRMSD();
        int init_obs();
        int init();
        //! Analyze Ramachandran angle differences from reference
        double evaluate();

        void avg_fill(int i);
        void avg_reset();
        void avg_write(Output &op);

        void reference_angle_list(std::vector<int> &resind,
                                  std::vector<double> &phvec,
                                  std::vector<double> &psvec);
        void reference_angle_file(std::string filename);
        int reference_angle_pdbfile(std::string filename);
        inline int num_used_res() const {return (int) refres.size();}

        //! Squared RC angle difference for the i'th residue in the whole system
        inline double rc_sq_diff(int i) const { return diff2res[i]; }

        void rangeEstimate(double &x0, double &x1);
    private:
        double diff2mean;
        std::vector<double> refphi,refpsi;
        std::vector<int> refres;
        std::vector<double> diff2res;
        std::valarray<unsigned long> dndt;
        Matrix<double> resprops;
    };
}

#endif
