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

#ifndef Ligand_HH
#define Ligand_HH
#include <deque>
#include "../Elements/GroupLib.hh"
#include "../Elements/Atom.hh"
#include "../Elements/Dipole.hh"
#include "../Elements/Node.hh"
#include "../Aux/ConnectionsMatrix.hh"
#include "../Aux/fileutils.hh"
#include "../Aux/AtomRecord.hh"
#include "../Aux/prf_xml.hh"

namespace prf
{
    //! An abstract class providing a base for amino acids and capping groups
    /**
    * The word "Ligand", as used inside PROFASI is an abstraction. It means a
    * collection of atoms with covalent connections. It has a name and a few
    * degrees of freedom, and it can be a part of a protein chain. Amino acids
    * are, in this sense, Ligands, and so are capping groups like Acetyl or NH2.
    * It makes no sense to define a bare, "Ligand" object, as it has few useful
    * properties. It will have atoms, but will not know how to place them in
    * space, for instance. So, the Ligand class does not have a public
    * constructor. One can only inherit from this class to create a more useful
    * object and then use the inheritted class in the program.
    *
    * \ingroup physical_entities
    */

    class Ligand
    {

    public:
        Ligand(const Ligand &);
        Ligand & operator=(const Ligand &);
        virtual ~ Ligand();
        virtual void Allocate();
        virtual void Initialize();
        inline OneLetterCode OLC() const {return mytype;}

        inline char CharCode() const {return Groups::mapOLC2Char(mytype);}

        inline std::string Name() const {return mygrp->CommonName();}

        inline std::string TLC() const {return mygrp->TLC();}

        inline void setName(OneLetterCode olc) {mytype = olc;}

        inline Atom & ATOM(int i) {return atm[i];}

        inline Atom & atom(int i) {return atm[i];}

        inline Atom & memberAtom(int i) {return atm[i];}

        inline Atom & first_atom() {return atm[0];}

        inline Atom & last_atom() {return atm[atm.size()-1];}

        inline int NumberOfAtoms() {return natms;}

        inline int numAtoms() {return natms;}

        inline int numHeavyAtoms() const {return nhvatm;}

        inline Ligand *RightConnection() {return rc;}

        inline Ligand *LeftConnection() {return lc;}

        inline void RightConnection(Ligand * glg) {rc = glg;}

        inline void LeftConnection(Ligand * glg) {lc = glg;}

        inline int UniqueId() const {return unid;}

        inline void UniqueId(int i) {unid = i;}

        inline int LocatedOn() const {return chloc;}

        inline void LocatedOn(int i) {chloc = i;}

        inline int atomOffset() const {return atOffset;}

        virtual void atomOffset(int i);
        inline int rotOffset() const {return rtOffset;}

        inline void rotOffset(int i) {rtOffset = i;}

        inline int n_rotDof() const {return nrtdof;}

        inline void n_rotDof(int i) {nrtdof = i;}

        virtual Node *node_for_dof(int i);

        //! Assign to (sidechain) dof i, and return the range of moved atoms
        /**
        * Assigns to the i'th degree of freedom of this Ligand and reconstructs.
        * Puts the unique ids of the first and the last changed atoms in a0 and a1.
        * An important difference from the older, obsolete, ROTDOF(...) is that the
        * index i is interpreted as the index within this Ligand, not the index
        * with respece to the whole system.
        */
        virtual int rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1);
        //! Assign to (sidechain) dof i, but don't recalculate any coordinates
        virtual int rotDof_assign(int il, double mgd);
        //! Value of the i'th (sidechain) dof
        virtual double get_rotDof(int il);
        //! Get the atoms defining the axis of the i'th sidechain dof
        virtual bool get_rotDofAxis(int i, Atom &a0, Atom &a1);

        //! Value of the i'th coordinate
        /**
        * For an EndGroup, this would mean the torsional degrees of freedom
        * in the group. For an AminoAcid, it would be the phi, psi, omega
        * angles, followed by the side chain degrees of freedom. Note that
        * the omega angles are not really degrees of freedom in PROFASI.
        * Neither is the phi angle for Proline or D-Proline. But for uniformity
        * of the backbone description, they are included in what are generically
        * called "coordinates". The specialized functions like BBDOF() in
        * Protein do not count omega or Proline phi angles in the list of
        * degrees of freedom. But Protein::get_coord(i) will return the omega
        * angle if i corresponds to it.
        */
        virtual double get_coord(int i);

        //! Set the i'th coordinate
        /**
        * Same as get_coord() but this function sets the value.
        */
        virtual void set_coord(int i, double x);
        //! Total number of coordinates
        /**
        * Number of degrees of freedom, including fixed backbone angles like
        * omega angles and proline phi angles.
        */
        virtual int n_coord() const;
        //! Get the value of the i'th true degree of freedom
        virtual double get_dof(int i);
        //! Set the value of the i'th true degree of freedom
        virtual void set_dof(int i, double x);
        //! Total number of true degrees of freedom
        virtual int n_dof() const;

        //! Assign degrees of freedom from an XML_Node object
        /**
        * This function tries to extract values for all the degrees of freedom
        * of the Ligand from an XML_Node, passed through its pointer. It checks
        * if the 3 letter code in the XML_Node matches with the 3 letter code of
        * the Ligand. If they match, dof info is extracted from the node and
        * assigned. If they don't and the flag typecheck is non-zero, an error
        * results. If typecheck is 0, the ligand type mismatch is ignored, and
        * the everything that can be assigned from the XML node is assigned.
        *
        * Note: At present XML tag "coordinates" and "dofs" are treated in the
        * same way in this function. The function can only read in coordinates,
        * which means, a node should have values for phi, psi, omega even if
        * the residue is a proline. The values passed as interpreted as
        * coordinates, not as degrees of freedom. It does not harm to pass
        * rubbish values for the angles which are not real degrees of freedom.
        * Nothing is ever assigned to torsion angles which are model constants
        * in ProFASi.
        */
        virtual int set_coord_xml(prf_xml::XML_Node *rs, int typecheck=1);

        inline bool isAA() const {return isaa;}

        inline bool isEG() const {return iseg;}

        inline bool isHydrophobic() const {return mygrp->isHydrophobic();}

        inline bool isPolar() const {return !isHydrophobic();}

        inline bool isCharged() const {return mygrp->isCharged();}

        inline int NumberOfDonors() const {return nd;}

        inline int NumberOfAcceptors() const {return na;}

        inline int SeqSerial() const {return seqser;}

        inline void SeqSerial(int i) {seqser = i;}

        inline int chargedAA_index() const {return chser;}

        inline void chargedAA_index(int i) {chser = i;}

        inline int n_donors() const {return nd;}

        inline int n_acceptors() const {return na;}

        virtual void Donor(int i, Dipole & dp);
        virtual void Acceptor(int i, Dipole & dp);
        virtual void ExportConnections(ConnectionsMatrix & aa);
        virtual void LocPairsatRTdof(int i,
                                     std::deque<std::pair<int,int> >&lcp);
        virtual std::string label_of(int i);
        virtual Atom * labeled_atom(std::string alabel);
        virtual Atom & at(std::string alabel);
        virtual void WritePDB(int &atindx, char ch_id, int aaindx,
                              FILE *fp);
        virtual void WritePDB2(int &atindx, char ch_id, int aaindx,
                               FILE *fp);
        void imprint(int iat,AtomDescriptor &ds);
        void WritePDBline(int &atindx, char ch_id, int aaindx,
                          int i, FILE *fp, int het=0);
        void Write_XML(FILE *fp);
        prf_xml::XML_Node * make_xml_node();
        int read_coordinates(std::list<AtomRecord>::iterator start,
                             std::list<AtomRecord>::iterator end,
                             std::vector<bool> &assignments);
        int calc_torsions(std::vector<bool> &specified);
        virtual void Reconstruct();
        virtual void BuildConnections();

        //OBSOLETE!! Will be removed from future versions. 2007-09-04
        virtual int ROTDOF(int il, double mgd, int &a0, int &a1);
        virtual int ROTDOFr(int il, double mgd);
        virtual double ADOF(int il);

    protected:
        //prevent construction of basic ligand
        explicit Ligand(OneLetterCode cod);
        Ligand *rc, *lc;
        int unid,chloc,atOffset,nrtdof,rtOffset,na,nd,natms,nhvatm;
        int seqser, chser;
        std::vector<Atom> atm;
        std::vector < Node * > node;
        bool isaa,iseg;
        OneLetterCode mytype;
        GroupProps * mygrp;
    };
}

#endif
