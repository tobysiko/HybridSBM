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

#ifndef HydrogenBond_HH
#define HydrogenBond_HH
#include <vector>
#include "../Elements/Population.hh"
#include "../Updates/Update.hh"
#include "../Elements/Dipole.hh"
#include "../Aux/Matrix.hh"
#include "Energy.hh"

namespace prf
{
//! HydrogenBond base class from which the HBMM and HBMS classes are derived
    /**
     * The HydrogenBond base class implements the basic functional form of the
     * hydrogen bond terms in PROFASI. The terms represented by HBMM and HBMS
     * make use of backbone NH and CO dipoles. It is convenient to store
     * information about them in one place.
     * \ingroup profasi_energies
     *
     */

    class HydrogenBond
    {
    public:
        HydrogenBond();
        virtual ~ HydrogenBond();
        virtual void set_population(Population * p);
        inline void EpsHB1(double x) { epshb1 = x; }

        inline void EpsHB2(double x) { epshb2 = x; }

    protected:
        //! Bare hydrogen bond function between two dipoles.
        /**
         hb_pair, in version 1.1 of profasi, calculates the interaction
        between a pair of dipoles. Dipoles now know which ligand and
        chain they are on. hp_pair makes use of this information to
        decide whether or not to use periodic boundary conditions. Since
        by contruction in profasi, no pair of atoms in a single protein
        can have a distance greater than one box length, periodic distance
        measure is not used if both dipoles are in the same chain. The
        ligand index is used to enforce the suppression of hydrogen
        bonds between neighbouring residues. The dipoles belonging to
        neighbouring peptide units do not participate in hydrogen bonding.
        For each dipole, there is therefore an excluded region in the list
        of dipoles. The region is given in relative ligand indices with
        respect to the ligand contributing the donor. If excl1 is -2 and
        excl2 is 0, it means, the contribution is 0 for acceptor ligand
        indices idon-2, idon-1 and idon.
        */
        double hb_pair(Dipole &don, Dipole &acc, int excl1,int excl2);
        double dhb_pair(Dipole &don, Dipole &acc, int excl1,int excl2,
                        std::valarray<double> &gx);
        static inline double sqr(double x) {return x * x;}

        double sighb, cuthb, sighb2, cut2, ahb, bhb, epshb1, epshb2;
        double POWA, POWB, halfpowa, halfpowb;
        double cdon, cacc, csacc, ntdonsupr,ntaccsupr,ctdonsupr,ctaccsupr;
        int ndon,nacc;
        std::vector<Dipole> donor,acceptor;
        std::vector<int> donstart, donend, accstart,accend;
    };

//! The Backbone-Backbone hydrogen bond term.
    /**
     * The hydrogen bonds represented by this class are formed between NH
     * and CO pairs from the backbones. A backup matrix is used to store
     * the contribution from the bond between donor i, and acceptor j.
     * The deltaE function is then optimized to recalculate only the changed
     * contributions. For instance when only a side chain torsional angle is
     * changed, nothing changes for the backbone hydrogen bonds. When a pivot
     * operation is carried out on one chain, backbone-backbone hydrogen
     * bonds in the moved part of that  chain, and their interaction with all
     * other chains are recalculated. But not the backbone-backbone
     * interactions which don't involve the moved parts. The newly calculated
     * terms are compared with the corresponding terms in the backup matrix
     * and the sum is the required energy change.
     * \ingroup profasi_energies
     *
     */

    class HBMM: public HydrogenBond, public Energy
    {
    public:
        HBMM();
        ~HBMM();
        void init();
        double evaluate();
        double gradientXYZ(std::valarray<double> &gx);
        double deltaE(Update *);
        void Accept(Update *);
        //! Inter chain hydrogen bonds (asymmetric)
        /**
        * Chain ich contributes NH groups, chain jch contributes CO groups
        */
        double InterChain(int ich, int jch);
        //! Value of the inter amino acid hydrogen bond potential
        /**
        * ich and jch are the chain indices. iaa and jaa are the
        * aminoacid indices inside those chains. This is an asymmetric
        * function. ich and iaa represent the aminoacid that contributes
        * the donor and jch and jaa represent the one that contributes
        * the acceptor. So, the value obtained when (ich,iaa) is replaced
        * by (jch,jaa) will in general be different.
        */
        double InterAA(int ich, int jch, int iaa, int jaa);
        //! Value of the inter residue hydrogen bond potential
        /**
        * Same as above, except that for each chain, the ligand indices
        * instead of the amino acid indices are used.
        */
        double InterLg(int ich, int jch, int ilg, int jlg);
        //! Inter amino acid hydrogen bond energy
        /**
        * For this function, iaa and jaa mean the ligand indices in the
        * whole system. That is, if there is more than one chain, residues
        * are numbered in an obvious way, by numbering residues from one
        * chain after the residues of another. This function is asymmetric
        * in iaa and jaa. The amino acid corresponding to iaa contributes
        * the donor, and jaa contributes the acceptor.
        */
        double InterLg(int iaa, int jaa);
        //! Inter amino acid hydrogen bond energy
        /**
        * As above but with amino acid indices relative to the whole system.
        */
        double InterAA(int iaa, int jaa);
        double interhb();
        void rangeEstimate(double &x1, double &x2);
    private:
        std::valarray<double> Mmm, Vmm;
        std::valarray<int> changed;
        int nchanges;
    };

//! The Backbone-Sidechain hydrogen bond term
    /**
     * This class represents the hydrogen bonds between the NH and CO dipoles of
     * the backbones with charged side chains of all proteins in the system. Only
     * amino acids D, E, K and R make side-chain backbone hydrogen bonds in the
     * model. The term is significantly weaker than the HBMM term. The
     * implementation here takes care of the optimization of this energy term. For
     * a side chain update on a charged side-chain, only the interaction of that
     * side chain with all the backbones is recalculated, and compared with a
     * backup matrix. When a backbone or rigid body update is made, interaction
     * of all side chains in other protein objects with that backbone, as well as
     * interactions of all side chains attached to that backbone with all other
     * backbones, are recalculated.
     * \ingroup profasi_energies
     *
     */

    class HBMS: public HydrogenBond, public Energy
    {
    public:
        HBMS();
        ~HBMS();
        void init();
        double evaluate();
        double gradientXYZ(std::valarray<double> &gx);
        double deltaE(Update *);
        void Accept(Update *);
        void check_contrib();
        void rangeEstimate(double &x1, double &x2);
    private:
        bool participates(prf::OneLetterCode cd);
        std::valarray<int> changed;
        int nchanges,nscdon,nscacc;
        std::vector<Dipole> scdonor,scacceptor;
        std::vector<int> scdonbeg,scaccbeg,scdonend,scaccend;
        std::valarray< double > Mms, Vms;
    };

}

#endif
