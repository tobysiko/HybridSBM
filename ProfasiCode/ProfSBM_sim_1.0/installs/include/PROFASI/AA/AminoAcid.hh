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

#ifndef AminoAcid_HH
#define AminoAcid_HH
#include "../Aux/Constants.hh"
#include "Ligand.hh"
#include "../Elements/Backbone.hh"
#include "../Aux/profasi_io.hh"
#include "../Aux/FindCoord.hh"
#include <string>

/**
* \defgroup aminoacids Amino Acids
* \ingroup physical_entities
*/

namespace prf
{

    typedef enum {LEV, DEX} ChiralityType;
//! Amino acid base class
    /**
     * This class implements most of the properties of the amino acids. The precise
     * geometry of individual amino acids is left to the derived classes.
     * \ingroup aminoacids
     */

    class AminoAcid: public Ligand
    {

    public:
        explicit AminoAcid(OneLetterCode typ);
        AminoAcid(const AminoAcid &) ;
        AminoAcid &operator=(const AminoAcid &);
        static void initCommon();
        virtual ~ AminoAcid();
        inline std::string SCALabel(int i) {return mygrp->label(icb+i);}

        //! Reference to backbone nitrogen
        inline Atom & Nitrogen() {return atm[0];}

        //! Reference to Calpha
        inline Atom & Calpha() {return atm[ibbca];}

        //! Reference to C' atom
        inline Atom & Cprime() {return atm[ibbc];}

        //! Reference to H alpha, or first H attached to C alpha
        inline Atom & Hca() {return atm[ibbca+1];}

        //! Reference to C beta
        inline Atom & Cbeta() {return atm[icb];}

        //! Reference to O attached to C'
        inline Atom & Oc() {return atm[ioc];}

        //! Reference to i'th side chain atom
        inline Atom & sidechain_atom(int i) {return atm[icb+i];}

        //! Reference to any atom in the amino acid through its PDB label
        /**
         * This is the slowest way to access an atom, and should be used with
         * care. It is harmless during initialization and for small analysis
         * programs, but disasterous to the performance if you decide to use
         * this, for instance, in an energy calculation. Note also that it
         * returns a pointer rather than a reference like the above functions.
         * This is necessary to address the cases when a labeled atom is asked
         * which does not exist, like CG2 on alanine. In such cases we want to
         * be able to return NULL.
         */
        Atom * labeled_atom(std::string alabel);
        //! Look up and return reference to Atom with a certain label
        /**
        * Warning: Upon failure to find the requested label, this function
        * can not return anything reasonable. So, if that happens,
        * the program exits.
        */
        Atom & at(std::string alabel);
        std::string label_of(int i);
        inline int numSideChainAtoms() const {return nsatm;}

        //! Side chain torsional DOF numbered from 0
        inline double Chi(int i) const {return node[i]->Phi();}

        //! Ramachandran Phi angle about the bond N-CA in the residue
        inline double Phi() const {return theBB->torsion_angle(BBloc);}

        //! Ramachandran Psi angle about the bond CA-C in the residue
        inline double Psi() const {return theBB->torsion_angle(BBloc+1);}

        //! Peptide bond angle between this and the next residue
        inline double Omega() const {return theBB->torsion_angle(BBloc+2);}

        //! Set side chain torsional DOF
        inline void Chi(int i, double x) { node[i]->AssignPhi(x); }

        //! Set Ramachandran Phi angle
        inline void Phi(double x) { theBB->torsion_angle(BBloc,x); }

        //! Set Ramachandran Psi angle
        inline void Psi(double x) { theBB->torsion_angle(BBloc+1,x); }

        //! Set omega angle between this and the next residue
        inline void Omega(double x) {
            if (cos(x)>=0) theBB->torsion_angle(BBloc+2,0);
            else theBB->torsion_angle(BBloc+2,UnivConstants::pi);
        }

        inline double Rotamer(int i) const {return node[i]->Phi();}

        inline double sc_torsional_angle(int i) const {return node[i]->Phi();}

        inline int numSideChainDOF() const {return nrtdof;}

        inline bool hasNTerminal() const {return ntrml;}

        inline bool hasCTerminal() const {return ctrml;}

        inline void hasNTerminal(bool trstat) {ntrml = trstat;}

        inline void hasCTerminal(bool trstat) {ctrml = trstat;}

        //returns true if the peptide bond between this residue and
        //the next is in cis-conformation
        inline bool isCis() const {return cis;}

        //set the succeeding peptide bond to cis
        inline void setCis(bool c = true) {cis = c;}

        virtual void Reconstruct();
        virtual int UpdateSideChainDof(int rdof);
        virtual void Initialize();
        inline ChiralityType Ca_Chirality() const {return chirality;}

        void Ca_Chirality(ChiralityType c) {chirality=c;}

        //! Calculate Ramachandran angles from current 3D coordinates.
        void calcPhiPsi(double &phv, double &psv) const;
        //! Do the current backbone angles fall in the helix region ?
        bool is_helical() const;
        //! Do the current backbone angles fall in the beta strand region ?
        bool is_strand() const;
        void Write();
        void WriteNodes();
        void WriteNodes(std::string &nds);
        //! Writes into a PDB file.
        /**
         * \c aaindx is the serial number to be used for the
        * amino acid. \c ch_id is the chain identifier.
        * \c atindx is the starting atom index.
         */
        void WritePDB(int &atindx, char ch_id, int aaindx, FILE * fp);
        //! Write PDB with heavy atoms first.
        void WritePDB2(int &atindx, char ch_id, int aaindx, FILE * fp);

        //! phi, psi, omega and side chain chi angles
        /**
        * For all amino acids, get_coord returns the angle phi for i=0, psi for
        * i=1 and the omega angle leading to the next residue for i=2. For
        * i=3,4,5... the side chain angles chi0, chi1 etc are returned. Note,
        * that omega is not really a degree of freedom in PROFASI, as of May
        * 2010. Similarly phi angle of Proline and D-Proline are fixed angles.
        * But for simplicity, they are returned here in their respective
        * places. get_dof(i) on the other hand, will skip those angles and
        * return the i'th true degree of freedom.
        */
        double get_coord(int i);
        //! Assign to i'th DOF in amino acid
        /**
        * Just like get_coord(), but this function assigns values.
        */
        void set_coord(int i, double x);
        //! 3+number of side chain angles
        /**
        * Backbone phi, psi and omega angles and side chain degrees of freedom
        * are counted.
        */
        int n_coord() const;
        //! True degrees of freedom in the residue
        /**
          * Returns phi,psi,chi0... for all non-proline residues and only psi
          * for proline.
        */
        double get_dof(int i);
        //! Assign to i'th true DOF in amino acid
        /**
        * Just like get_dof(), but this function assigns values.
        */
        void set_dof(int i, double x);
        //! Number of true degrees of freedom in the residue
        /**
        * Backbone phi (non-proline), psi and side chain degrees of freedom
        * are counted.
        */
        int n_dof() const;

//The functions that follow are only for internal use inside profasi. They
//should not be used in application programs.
        inline int locIndCa() const {return ibbca;}

        inline int locIndCb() const {return icb;}

        int rotDof_assign_and_reconstruct(int il, double mgd, int &a0, int &a1);
        int rotDof_assign(int il, double mgd);
        double get_rotDof(int il);

        int ROTDOF(int il, double mgd, int &a0, int &a1);
        int ROTDOFr(int il, double mgd);
        double ADOF(int il);
        std::string USCAtomDescr(int iscatm);
        void Donor(int i, Dipole & dp);
        void Acceptor(int i, Dipole & dp);
        void ExportConnections(ConnectionsMatrix & aa);
        void BuildConnections();
        void LocPairsatRTdof(int i, std::deque < std::pair < int,int > >&lcp);
        inline void AttachBB(Backbone *gbb, int attachpt) {
            theBB=gbb;
            BBloc=3*attachpt;
        }

        void Allocate();
        //! Reassign atom offsets to the nodes
        /**
        * Nodes maintain their own atom lists. If atoms of an amino acid get
        * an offset, the nodes are no longer connected to the right atoms.
        * This function restores the connectivity.
        */
        void nodes_reconnect();
        static const double theta[3];
        static const double b[3];
        static const double bNH, bCaHa, bCaCb, bCO;
        static const double phCb, thCb, phHa, phHa2, thHa,DPRphHa,DPRphCb;
        static FindCoord locate_Hn, locate_Ha, locate_Ha2,
        locate_Oc, locate_Cb, locate_O2c, locate_H1n,
        locate_H2n, locate_H3n,locate_dHa,locate_dCb;
        static bool initialized;

    protected:
// void AddGroup(int &i, ConstituentGroup gp);
        inline void Add(int &i, AtomKind att) {atm[i++].Species(att);}

        inline void charCode(char ch) {mytype = Groups::mapChar2OLC(ch);}

        bool ntrml, ctrml, cis;
        int nsatm, ibbca, ibbc, ioc , icb;
/*        Vector3 *bcn1, *bnca, *bcac, *bcn2;*/
        Backbone *theBB;
        int nnodes, BBloc;
        ChiralityType chirality;
    };
}

#endif
