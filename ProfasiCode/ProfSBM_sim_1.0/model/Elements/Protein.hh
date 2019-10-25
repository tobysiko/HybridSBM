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

#ifndef Protein_HH
#define Protein_HH
#include <vector>
#include <set>
#include "../Aux/RandomNumberBase.hh"
#include "NaturalAminoAcids.hh"
#include "Backbone.hh"
#include "KnownEndGroups.hh"
#include "../Aux/ConnectionsMatrix.hh"

namespace prf
{
    //! Representation for a poly peptide chain
    /**
      The Protein class in ProFASi represents a string of amino acids
      joined by peptide bonds. For the purposes of this package, the word
      protein is synonymous with a peptide or a poly-peptide chain. The
      molecule may contain non-amino acid capping groups (See EndGroup)
      at one or both ends of the chain.

      Geometrically, the molecule is represented as a sequence of short trees
      atoms connected to a Backbone object. The Protein class provides
      information about the constituent parts and the degrees of freedom
      of the protein as well as providing functions for their construction and
      manipulation.
      \ingroup physical_entities
    */

    class Protein
    {
    public:
        //! Default constructor
        Protein();
        //! Constructor taking an amino acid sequence
        Protein(std::string aaseq);
        //! Copy constructor
        Protein(const Protein &);
        //! Assignment operator
        Protein & operator=(const Protein &);
        ~Protein();
        //! Translate molecule
        void translate(const Vector3 &trv)
        {
            AtomCoordinates::BlockTranslate(trv,atbg,atnd);
        }

        //! Rotate molecule about an axis passing through the center of mass.
        void rotate(double angl, const Vector3 &axs)
        {
            AtomCoordinates::BlockRotate(angl,atbg,atnd,CenterOfMass(),axs);
        }
        //! Recalculate cartessian coordinates from internal
        void reconstruct();
        //! Reconstruct given a starting point and a direction
        /**
         * reconstruct starting from amino acid iaa, to the C terminus if idir
         * is 0, to the N terminus if idir is 1. The rest of the chain is left
         * unchanged.
         */
        void reconstruct(int idir, int iaa);
        //! Bring all coordinates to their appropriate periodic intervals
        void EnforceBC();

        //@{
        /**
          @name Initialisation related functions
          */
        //! Clears all arrays and resets scalars to defaults
        /**
        * Deletes all pointers owned by the protein in pointer arrays and then
        * clears the arrays. Stored numbers such as the number of hydrophobic
        * amino acids are zeroed. Essentially it cleans all information stored
        * in a protein object, so that a new allocation process can proceed
        * without trouble.
        */
        void clear();
        //! Increment amino acid sequence with one residue at the C-terminal
        /**
        * This function appends an amino acid of a given type at the end of
        * the amino acid chain, on the C-terminal side. If there are end groups
        * attached to the protein already, this function makes sure to shift
        * the connection of the C-terminal  end group appropriately. Sometimes
        * it adjustd the number of atoms in the last amino acid before this>
        * addition (if it had an OXT, it has to lose it). The information
        * about the degrees of freedom and other stuff such as the arrays of
        * hydrophobic and charged residues are updated.
        *
        * This function changes the number of atoms in the protein. This will
        * in general require reassignment of the UniqueIds of all atoms in the
        * population, and a resizing of the atom coordinate storage. Atom
        * coordinates are not owned by proteins or residues or ligands or even
        * atoms!
        */
        int add_aminoacid(prf::OneLetterCode cod);
        //! Add N-terminal capping group
        int add_NTLigand(prf::OneLetterCode cod);
        //! Add C-terminal capping group
        int add_CTLigand(prf::OneLetterCode cod);
        //! Set a new sequence and allocate memory accordingly.
        /**
        * The end result of this function should be a properly initialized
        * Protein  object with a sequence as provided as an argument. If the
        * required sequence differs from the actual sequence, the function
        * "clear()" is called, followed by a sereis of calls to add_aminoacid
        * and add_..._endgroup as required.
        *
        */
        int set_sequence(const std::vector<prf::OneLetterCode> &gseq);
        void connect_residues();
        void init_atoms();
        void create_backbone();

        //! Initialize DOF arrays
        void init_dof();
        //! All internal coordinates to zero
        void trivial_init();
        //! Initialize to stretched out state.
        /**
         * Backbone angles are put to (-pi,pi) throughout (with the exception
         * of proline for which the backbone phi angle is fixed) and side
         * chain angles are put to 0.
         */
        void stretched_init();
        //! all internal coordinates to random values
        inline void random_init(RandomNumberBase *rangen) {randomize(rangen);}
        void randomize(RandomNumberBase *rangen);
        //! ideal alpha helix for backbone, zero for rotamers
        void helical_init();
        //! Adopt the sequence given in XML node.
        /**
        * The protein checks if there is a sensible sequence given in the
        * XML node. Checks if the sequence is already the same as its own.
        * If it discovers a new sensible sequence in the XML Node, it
        * discards its current sequence and allocates amino acids and capping
        * groups as per the new sequence. The return value is 1 if something
        * changed in the protein. This function does not recover or assign
        * any degrees of freedom from the XML node. For that, use assign_dof().
        */
        int metamorphose(prf_xml::XML_Node * px);
        //! Set the peptide bond between iaa and iaa+1 to cis
        inline void setCis(int iaa) {cisres.push_back(iaa);}

        //! Whether or not to use charged chain ends
        /**
        * First, we note that whatever you set here is useless if there is
        * a real end group specified. Second, if there are no end groups,
        * by default, both ends are charged. So, if that is what you want,
        * don't do anything. If you are not using capping groups and yet
        * want neutral chain end/s, use this function to set which ends
        * are uncharged, by passing a "false" for that end.
        */
        inline void charged_end(bool ntcg, bool ctcg) {charged_NT=ntcg;charged_CT=ctcg;}
        //@}

        //@{
        /**
          @name Accessing constituent properties
          */

        //! Each protein in the population has a unique integer id
        inline int Id() const {return unid;}
        //! Assign integer id
        inline void Id(int i) {unid=i;}
        //! Retrieve sequence string
        std::string Sequence() const;
        //! Retreive OneLetterCode of the i'th residue
        inline OneLetterCode Seq(int i) const {return seq[i];}
        //! Number of amino acids
        inline int numAminoAcids() const {return naa;}
        //! Number of a.a.+ capping groups
        inline int numLigands() const {return lg.size();}
        //! Get N terminal capping group
        inline EndGroup * NtermLigand() {return NtLigand;}
        //! Get C terminal capping group
        inline EndGroup * CtermLigand() {return CtLigand;}
        //! Pointer to i'th amino acid
        inline AminoAcid * memberAA(int i) {return AAv[i];}
        //! Pointer to i'th amino acid
        inline AminoAcid * AA(int i) {return AAv[i];}

        //! First ligand in chain
        inline Ligand * first_ligand() {return lg.front();}
        //! Last ligand in chain
        inline Ligand * last_ligand() {return lg.back();}
        //! Same as first_ligand
        inline int begin_ligand() const {return 0;}
        //! One past the last ligand
        inline int end_ligand() const {return lg.size();}
        //! First amino acid in chain
        inline AminoAcid * first_AA() {return AAv.front();}
        //! Last amino acid in chain
        inline AminoAcid * last_AA() {return AAv.back();}

        //! i'th ligand, counting a.a. and capping groups
        inline Ligand * memberLigand(int i) {return lg[i];}

        //! Pointer to the backbone
        inline Backbone *backbone() {return &bb;}
        //! i'th backbone atom
        inline Atom & backbone_atom(int i) {return bb.atom(i);}
        //! i'th backbone atom
        inline const Atom & backbone_atom(int i) const {return bb.atom(i);}
        //! Number of atoms
        inline int numberOfAtoms() const {return natms;}
        //! First atom in terms of unique ids
        inline int begin_atom() const {return atbg;}
        //! One past the last element
        inline int end_atom() const {return atnd;}
        //! Same as begin
        inline int first_atom() const {return atbg;}
        //! Last atom
        /**
         * end_atom is like the end() function of standard library containers,
         * one past the last element, to be used in loops with a "<" condition
         * for termination. last_atom is really the last atom. Use it with
         * "<=" terminating condition in loops. The integer return values are
         * UniqueIds
         */
        inline int last_atom() const {return atnd-1;}
        //@}

        //@{
        /**
          @name Miscellaneous properties of the molecule
          */

        double helix_fraction();
        double strand_fraction();

        //! End-to-end length
        double EndToEndDistance() const {return bb.EndToEndLength();}
        //! End-to-end vector
        Vector3 EndToEndVector() const {return bb.Axis();}
        //! Center of mass
        Vector3 CenterOfMass() const {
            return AtomCoordinates::CenterOfMass(atbg,atnd);
        }
        //! Squared backbone rd. of gyration
        double BackboneRg2() const {return bb.Radius2OfGyration();}

        //@}

        //@{
        /**
          @name Accessing and manipulating degrees of freedom
          */

        //! Access Ramachandran phi angle of the i'th amino acid
        inline double RamachandranPhi(int i)const {return bb.torsional_angle(3*i);}
        //! Access Ramachandran psi angle of the i'th amino acid
        inline double RamachandranPsi(int i)const {return bb.torsional_angle(3*i+1);}

        //! Assign to the Ramachandran phi angle of the i'th residue
        inline void RamachandranPhi(int i,double xv) {bb.torsional_angle(3*i,xv);}
        //! Assign to the i'th Ramachandran psi angle
        inline void RamachandranPsi(int i,double xv) {bb.torsional_angle(3*i+1,xv);}

        //! Number of (slightly redundant) degrees of freedom
        /**
          The return value is exactly 3 more than the true number of degrees
          of freedom for a protein in ProFASi. It is the sum of the number of
          backbone DOFs, side chain DOFs and 9 numbers to specify the overall
          position and orientation of the molecule. 6 numbers should be enough
          for the rigid body coordinates, but in ProFASi, we use 9, and that
          is the redundancy. The 9 numbers used are the Cartesian coordinates
          of the N, C_alpha and C' atoms of the first amino acid. Specifying
          these along with the torsional coordinates fixes the shape, position
          and orientation of the protein in space. Specifying the center of
          mass position and three Euler angles would suffice, but that involves
          more operations whenever these degrees of freedom need to be
          translated into Cartesian coordinates for all atoms.
          */
        int n_dof() const {return 9+numBBdof()+numRTdof();}
        //! Number of generalised coordinates
        /**
          This is similar to the n_dof function above, but here even the
          backbone torsion angles kept fixed in ProFASi are counted. For
          instance, the proline phi angle and the omega angles are kept fixed
          in ProFASi. So, for the sequence GPA, n_dof would return 9+5+1=15,
          where as this function would return, 9+9+1=19.
          */
        int n_coords() const {return 9+3*naa+numRTdof();}

        //! Size in bytes of the data block in ProFASi binary conf file
        inline int ConfSize() const {
            return (sizeof(double))*(9+numBBdof()+numRTdof())/sizeof(char);
        }
        //! Get the value of the i'th true degree of freedom
        double get_dof(int i);
        //! Set the value of the i'th true degree of freedom
        void set_dof(int i, double vl);
        //! Get the value of the i'th coordinate
        double get_coord(int i);
        //! Set the value of the i'th coordinate
        void set_coord(int i, double vl);

        //! Retrieve all dof in an array
        void get_dof(std::vector<double> &vdof) const;
        //! Retrieve all coordinates in an array
        void get_coord(std::vector<double> &vcrd) const;

        //! Set all dof from an array
        inline void set_dof(std::vector<double> &vdof)
        {set_dof(vdof.begin(),vdof.end());}
        inline void set_coord(std::vector<double> &vcrd)
        {set_coord(vcrd.begin(),vcrd.end());}

        void set_dof(std::vector<double>::iterator c1,
                     std::vector<double>::iterator c2);
        void set_coord(std::vector<double>::iterator c1,
                       std::vector<double>::iterator c2);

        //! Total number of backbone DOF
        inline int numBBdof() const {return bb.numDOF();}
        //! Total number of side chain DOF
        inline int numRTdof() const {return nrtval;}
        //! Access i'th backbone DOF
        inline double BBdof(int i) const {return bb.DOF(i);}
        //! Access i'th side chain DOF
        inline double RTdof(int i) const {return lg[rtloc[i]]->ADOF(i);}
        //! Assign to i'th backbone DOF
        inline void BBdof(int i, double val) {bb.DOF(i,val);}
        //! Assign to i'th side chain DOF
        inline int RTdof(int i, double val) {return lg[rtloc[i]]->ROTDOFr(i,val);}

        //! The residue which has the i'th side chain DOF
        inline Ligand *residue_with_rt_dof(int i) { return lg[rtloc[i]]; }
        //! The residue which has the i'th backbone DOF
        inline Ligand *residue_with_bb_dof(int i) { return AAv[bb.DOFloc(i)]; }

        //@}
        //@{
        /**
          @name Writing structure in various formats
          */
        //! Write down the molecule in plain text, including coordinates etc
        void Write();
        //! Write PDB with the normal order of atoms
        void WritePDB(int &istatm,FILE *fp, char ch_id, int &rsindx);
        //! Write pdb with heavy atoms first for each amino acid.
        void WritePDB2(int &istatm,FILE *fp, char ch_id, int &rsindx);
        //! Write an XML Node in the ProFASi XML structure format
        void Write_XML(FILE *op);
        //! Retrieve XML_Node object with a map of the molecule
        prf_xml::XML_Node * make_xml_node();

        //! Signature of the data to be written in a binary conf file
        std::string ConfSignature();
        //! Write current state in compressed binary format
        void WriteConf(FILE *fp);
        //! Write bare minimum info about current state in plain text
        void WriteConf_text(FILE *fp);
        //@}
        //@{
        /**
          @name Reading in the structure
          */
        //! Read in a block from a binary conf file
        void ReadConf(FILE *fp);
        //! Read in coordinates from an open text file
        void ReadConf_text(FILE *fp);
        //! Copy Cartesian coordinates from a list of AtomRecord objects
        /**
         * A typical use of this function would be, to make a selection
         * of some part of a PDB file with a PDBReader object, and then
         * pass the generated list of AtomRecords to the Protein object
         * to give it the same structure. It does not try to adjust the
         * geometry to the model geometry or fill in unspecified coordinates.
         * This function can be called several times to assign coordinates
         * to different parts of the chain.
         *
         * Think of it as a list copy operation. The protein is like a
         * list of atoms, and it imports entries from another list from
         * within limits specified by the start and end iterators. The
         * residue identifiers in the AtomRecord list are used only to
         * identify blocks meant for one residue, i.e., their numerical
         * values are ignored. The copying commences at a given location
         * along the chain \c at_ligand, and successive blocks in
         * the AtomRecord list are assigned to successive residues. The
         * \c assignments is used to record which atoms were actually
         * assigned to. It must be pre-initialized with a size equal to
         * the number of atoms in the entire population. For each atom
         * receiving coordinates through this function, the entry
         * corresponding to its unique id is set to "true". The entries
         * for the other atoms are not touched.
         * */
        int importXYZ(std::list<AtomRecord>::iterator start,
                            std::list<AtomRecord>::iterator end,
                            std::vector<bool> &assignments, int at_ligand=0);
        //! Calculate all internal DOF from Cartesian coordinates
        /**
          * This calculates the torsional degrees of freedom of the protein
          * from the Cartesian coordinates of the atoms.
          */
        int calc_torsions(std::vector<bool> &specified);

        //! Guess the coordinates of a few atoms based on other atoms
        /**
         In full generality, this function should try to guess the coordinates
         of some atoms if some of the atoms have been assigned proper coordinates.
         The input argument, a vector<bool> contains a list saying which
        coordinates are specified, which are not. The indices of the vector
        are supposed to be UniqueIds of the atoms. So, the size of the vector
        is the total number of atoms in the population. If, for instance,
        the position of all non-hydrogen atoms are specified,  simple
        triangulation can be used to fill in the hydrogen coordinates as
        best as possible. When parts of a side chain are missing but the
        degrees of freedom can still be determined it is possible to
        calculate the position of the missing atoms. If the degrees of
        freedom can not be determined, they can be set to their nominal
        values so that the side chain would not look ridiculous if it
        was alone. But if  the backbone atoms are missing, currently
        there is no remedy. So, this is intended to be used when a
        structure is imported from a PDB file or from  a reduced model,
        where the hydrogen coordinates are not specified, and the side
        chain atoms may not all be given. A little energy minimisation
        after this should give a workable model approximation to the
        imported structure.
        */
        int guess_missing_coordinates(std::vector<bool> & specified);

        //! Guess appropriate coordinates for a section of the backbone
        /**
          This assumes that X,Y,Z coordinates were read for atoms from
          somewhere (e.g., PDB file), but there were a few residues where
          even the backbone atoms were not specified. This function
          should close the chain by assigning coordinates to the
          unspecified part of the chain such that the model constraints
          are satisfied and the coordinates of the atoms which have been
          specified are not touched.
          @note At present this is a dummy. No algorithm to fix broken
          backbone is implemented in reality.
          */
        int fix_broken_backbone(std::vector<bool> &specified);

        //! Read in dof info from an XML node
        /**
        * Checks sequence against sequence stored in the XML node. If
        * they match, it extracts the child nodes with tag "group" and
        * assigns them to the appropriate residues. The nodes must have
        * an attribute "index", which stores the serial number for that
        * group in the protein sequence, counting from 0. It is possible
        * to have only a few "group" nodes in the XML file. One may be
        * interested in assigning coordinates to only a small block in the
        * middle, for instace. The entire sequence string is available in
        * the protein node of the XML file. The ligand indices are
        * interpreted with respect to this entire protein sequence.
        *
        * The optional flag "mismatch_strategy" tells the function what to do
        * in case the sequence in the XML node does not match the sequence
        * of the protein. Value 0 means, quit with an error message. This
        * is the default. A value of 1 means, assign whenever possible.
        */
        int assign_dof(prf_xml::XML_Node *px, int mismatch_strategy=0);
        //@}
    public:
//The following functions are used internally in profasi, but are not meant
//to be directly used by the end user.
        inline int BBAxisAtoms(Atom &at1, Atom &at2,int il)
        {return bb.AxisAtoms(at1,at2,il);}

        inline double RotDOF(int il) {return lg[rtloc[il]]->ADOF(il);}

        inline int RotDOF(int il,double mgd,int &at0, int &at1)
        {return lg[rtloc[il]]->ROTDOF(il,mgd,at0,at1);}

        inline int RotDOFr(int il,double mgd) {return lg[rtloc[il]]->ROTDOFr(il,mgd);}

        void Allocate(int ppid, int atom_offset,std::vector<OneLetterCode> fulseq);
        void Allocate(int ppid, int atom_offset,std::string sq);
        void LocPairsBBdof(int i, std::deque<std::pair<int,int> > & lcp);
        void LocPairsRTdof(int i, std::deque<std::pair<int,int> > & lcp);
        void BuildConnections();
        void ExportConnections(ConnectionsMatrix &connected);
        void setAtomOffset(int i);
        inline int AtomOffset() const {return atoffset;}
        void setUniqueId(int i);
        inline bool internally_consistent() { return consi; }
        bool isSameTypeAs(const Protein &p2);


    private:
        //UniqueIds of N,Ca and C for amino acid i.
        int iN(int i) {return AAv[i]->Nitrogen().UniqueId();}
        int iCa(int i) {return AAv[i]->Calpha().UniqueId();}
        int iC(int i) {return AAv[i]->Cprime().UniqueId();}
    private:
        std::vector<AminoAcid *> AAv;
        std::vector<Ligand *> lg;
        Backbone bb;
        EndGroup * NtLigand,*CtLigand;
        int unid;
        std::vector<OneLetterCode> seq;
        std::vector<int> rtloc,cisres;
        int atbg,atnd,lgbg,lgnd;
        int nrtval,natms,atoffset,naa;
        bool charged_NT, charged_CT, consi;
        void ExportConnectionssc(ConnectionsMatrix &connected);
        bool good_sequence(std::vector<OneLetterCode> sq);
    };

    typedef Protein Peptide;
    typedef Protein PolypeptideChain;
}
#endif
