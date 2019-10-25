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

#ifndef Backbone_HH
#define Backbone_HH
#include <vector>
#include <cmath>
#include "../Aux/Constants.hh"
#include "../Aux/FindCoord.hh"
#include "Atom.hh"
#include "../Aux/RandomNumberBase.hh"
#include "../Aux/Shape.hh"

namespace prf
{
    //Forward declaration of AminoAcid
    class AminoAcid;
    //! Representation for the protein backbone
    /**
     * \ingroup building_blocks
     * The Protein backbone is recognized as an important concept, and is
     * implemented in its own class. The atoms on the backbone of course
     * belong to different amino acids, and in principle it is possible to
     * represent the protein as just a chain of amino acids, with each
     * residue being responsible for constructing both the backbone and side
     * chain atoms. Practically however, it is simpler to construct the chain
     * make conformational updates on it ... if the backbone is constructed
     * together. Besides, any slight numerical precision errors in the
     * coordinates of the backbone atoms propagates much further than an error
     * of similar size in the placement of a side chain atom. So, in PROFASI
     * a very conservative, but safe reconstruction method is used to place
     * the backbone atoms in space.
     */

    class Backbone
    {
    public:
        Backbone(int naaa);   //!< Backbone for protein of naaa residues
        Backbone(const Backbone &);   //!< create a clone.
        ~Backbone();
        //! Clear all data for possible new initialisation
        void clear();
        //! Initialize common properties (eg. bond lengths..)
        void initCommon();
        void numAminoAcids(int naaa);   //!< Change number of residues to naaa
        //! Retrieve number of amino acids
        inline int numAminoAcids() const {return naa;}

        //! Number of backbone atoms
        inline int numAtoms() {return 3*naa;}

        //! Retrieve number of dof on backbone (taking into account proline..)
        inline int numDOF() const {return ndof;}

        //! Reference to i'th backbone atom
        inline Atom & memberAtom(int i) {return bbatom[i];}

        //! Reference to i'th backbone atom
        inline Atom & atom(int i) {return bbatom[i];}

        inline const Atom & atom(int i) const {return bbatom[i];}

        inline Vector3 *Bond(int i) {return &the_bond[i];}

        inline Vector3 bond(int i) {return the_bond[i];}

        //! N-terminus to C-terminus vector, from N to C
        inline Vector3 Axis() const
        {return bbatom[bbatom.size()-1].Pos()-bbatom[0].Pos();}

        //! End to end distance
        inline double EndToEndLength() const {return Axis().mag();}

        //! Backbone radius of gyration
        double Radius2OfGyration() const;
        //! Torsional angle along backbone, including the fixed ones
        inline double torsional_angle(int i) const {return phi[i];}
        inline double torsion_angle(int i) const {return phi[i];}

        //! Torsional angle along backbone, including the fixed ones
        inline void torsional_angle(int i,double xv) {if (!frzn[i]) phi[i]=xv;}
        inline void torsion_angle(int i,double xv) {if (!frzn[i]) phi[i]=xv;}

        //! Blind increment of a torsional angle, read details...
        /**
         * torsional angles means any torsional angle, even the fixed ones,
         * like the proline phi angle. DOF below means only those which
         * constitute degrees of freedoms for the molecule
         */
        inline void incr_tors_angle(int i, double val) {phi[i]+=val;}

        //! Access i'th degree of freedom along the backbone
        /**
         * The DOF functions take care to skip the fixed torsional angles
         * like the peptide bond angles (180 degrees) and proline phi angles
         */
        inline double DOF(int i) const {return phi[bbdof[i]];}

        //! Assign to the i'th dof
        inline void DOF(int i, double val) {phi[bbdof[i]]=val;}

        //! picks a DOF given a double between 0 and 1. Don't use.
        inline int DOFno(double lx) {return (int)(lx*ndof);}

        //! returns the index of the amino-acid that contains the i'th DOF
        inline int DOFloc(int i) const {return bbdof[i%ndof]/3;}

        //! returns the index in the list of torsional angles of the ith DOF.
        inline int DOFid(int i) const {return bbdof[i%ndof];}

        //! returns the index in the DOF list of ith torsional angle.
        inline int PhiToDOF(int i) const {return phitodof[i];}

        //! increment DOF i
        inline int incrDOFno(int i,double xval) {phi[bbdof[i]]+=xval;return 0;}

        //! increment a DOF corresponding to fraction between 0 and 1
        int incrDOF(double lx,double xval);
        //! Pick amino acid index based on random number xx. Don't use.
        int LocateAA(double xx);
        //! Find the two backbones connected by the bond with the i'th dof
        /** picks the two atoms connected by the bond which corresponds to the
         * DOF iloc. the return value is 1 if iloc is in the first half of the
         * backbone, else 0
         */
        int AxisAtoms(Atom &a0,Atom &a1,int iloc);
        inline int AxisAtoms(double xlocn,Atom &a0,Atom &a1)
        {return AxisAtoms(a0,a1,(int)(ndof*xlocn));}

        //! Register or connect an amino acid to the backbone
        /**
         * Backbone is created knowing only with the total number of amino
         * acids. This function is used to assign appropriate atoms to it, and
         * also to create a map between torsional angles and DOFs. For instance
         * if the amino acid is proline, the phi angle is not a DOF.
         */
        void registerAA(int i, AminoAcid *aac);
        //! Set atom offset to a given value
        /**
        * When a chain is created, it creates its own backbone. Start of the
        * backbone would then coincide with the N of the first amino acid. If
        * the chain is subsequently given an offset (if it is not the first
        * chain), the atoms of the backbone have to follow. This is the purpose
        * of this function. It is assumed that the relationship between the
        * UniqueIds of the atoms in one back bone do not change. You pass
        * the UniqueId of the first N after the offset. The other ids are
        * shifted appropriately.
        *
        * This was not necessary before, as the backbone was constructed once,
        * and the atoms in the protein already had proper UniqueIds when
        * the backbone was constructed. It is no longer true. A protein can
        * be deleted and reconstructed. In such case, the UniqueIds have to be
        * reassigned throughout the population after all protein deletion and
        * additions have been carried out.
        *
        */
        void set_atom_offset(int new_n0_uid);
        //! Randomize state of backbone using random number generator passed
        void randomize(RandomNumberBase *rangen);
        //! Set backbone state to a perfect alpha helix
        void SetToAlphaHelix();
        //! Set backbone state to a beta strand
        void SetToBetaStrand();
        //! calculate coordinates of backbone atoms in the N-C direction
        void forwardReconstruct(int strt, int iend);
        //! calculate coordinates of backbone atoms in the C-N direction
        void reverseReconstruct(int strt, int iend);
        //@{
        /**
         * @name Coordinates for the first and last three atoms
         * Assign coordinates to the first or last three atoms of the backbone.
         * This requires special treatment. When the backbone is reconstructed
         * from N-terminal to C-terminal, the positions of the first three
         * atoms are not fixed by the internal degrees of freedom of the
         * backbone. Similarly for the last three atoms. These functions
         * assign values near the origin if no previous value was assigned. If
         * the atoms are found to have valid positions, the functions will try
         * to keep them there while fixing round off errors in bond lengths
         * and the theta angle
         */
        void reconstStart();
        void reconstLast();
        //@}
        //! Bring all angles to the interval 0..2pi
        void FixAngles();
        //! Sanity check
        void CheckProperties();
        //! Find the i'th phi angle using the cartessian coordinates as input
        double find_phi(int i);
        //! Calculate and assign torsion angles from the Cartessian coordinates
        int calc_torsions(std::vector<bool> &specified);
        //! Find the i'th theta angle using the cartessian coordinates as input
        double find_theta(int i);
        //! Reconstruct bond vectors based on current coordinates of backbone atoms
        void reconst_bond_vectors();
        //! Reconstruct all backbone bonds to the right of backbone atom i1
        void freconst_bond_vectors(unsigned int i1);
        //! Reconstruct all backbone bonds to the left of backbone atom i1
        void breconst_bond_vectors(unsigned int i1);
        //! Export coordinates to a Shape object for use with RMSD
        void ExportCrds(Shape &sh) const;
        //! Export only part of the backbone for RMSD calculation
        void ExportCrds(Shape &sh, int iaast, int iaand) const;
    private:
        std::vector<Atom> bbatom;
        std::vector<Vector3> the_bond;
        std::vector<double> phi;
        std::vector<int> bbdof,phitodof;
        std::vector<bool> frzn;
        int naa,ndof,nmid;
        double b[3],theta[3];
        FindCoord bbb[3],bbf[3];
    };
}

#endif
