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

#ifndef PhiPsiChi_HH
#define PhiPsiChi_HH

namespace prf
{
    //! A single Ramachandran Phi angle as an observable
    /**
     * Useful when behaviour of individual degrees of freedom are interesting.
     * The only occasion when we used them was when we made maps for the
     * Ramachandran angles from our model for all amino acids.
     * \ingroup profasi_observables
     */

    class RC_phi : public Observable
    {
    public:
        RC_phi() : iaa(0),ch(NULL) {Name("Phi");}

        ~RC_phi() {}

        //! attach to some protein
        inline void set_protein(Protein *pr) {ch=pr;}

        //! attach to one amino acid
        inline void set_aa(int i) {iaa=i;}

        double evaluate() {
            return ch->RamachandranPhi(iaa);
        }

        void rangeEstimate(double &xmin,double &xmax) {
            xmin=-pi;xmax=pi;
        }

    private:
        int iaa;
        Protein * ch;
    };

    //! A single Ramachandran Psi angle as an observable
    /**
     * Useful when behaviour of individual degrees of freedom are interesting.
     * The only occasion when we used them was when we made maps for the
     * Ramachandran angles from our model for all amino acids.
     * \ingroup profasi_observables
     */

    class RC_psi : public Observable
    {
    public:
        RC_psi() : iaa(0),ch(NULL) {Name("Psi");}

        ~RC_psi() {}

        //! attach to a protein
        inline void set_protein(Protein *pr) {ch=pr;}

        //! attach to one amino acid
        inline void set_aa(int i) {iaa=i;}

        double evaluate() {
            return ch->RamachandranPsi(iaa);
        }

        void rangeEstimate(double &xmin,double &xmax) {
            xmin=-pi;xmax=pi;
        }

    private:
        int iaa;
        Protein * ch;
    };

    //! A single side chain Chi angle as an observable
    /**
     * Useful when behaviour of individual degrees of freedom are interesting.
     * The only occasion when we used them was when we made maps for the
     * Ramachandran angles from our model for all amino acids.
     * \ingroup profasi_observables
     */

    class Chi : public Observable
    {
    public:
        Chi() {Name("Chi");}

        ~Chi() {}

        //! attach to one protein
        inline void set_protein(Protein *pr) {ch=pr;}

        //! attach to one amino acid
        inline void set_aa(int i) {iaa=i;}

        //! specify which side chain angle
        inline void bind_to_dof(int i) {ichi=i;}

        double evaluate() {
            return ch->AA(iaa)->Chi(ichi);
        }

        void rangeEstimate(double &xmin,double &xmax) {
            xmin=0;xmax=2*pi;
        }

    private:
        int iaa,ichi;
        Protein * ch;
    };
}

#endif
