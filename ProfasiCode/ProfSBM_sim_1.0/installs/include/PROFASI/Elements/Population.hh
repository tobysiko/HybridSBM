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

#ifndef POPULATION_HH
#define POPULATION_HH
#include "Protein.hh"
#include "../Aux/Shape.hh"
#include "../Aux/AtomRecord.hh"
#include <deque>
#include <map>
#include "PopBase.hh"
#include "DOF_Index.hh"
#include "../Aux/prf_xml.hh"

namespace prf
{
    //! A population of proteins
    /**
     * Population is a collection of one or more Proteins of one or more kinds.
     * This is the class the conformational updates work on. This is the class
     * the energy terms calculate energies for. Population provides a convenient
     * interface to talk about "the system as a whole".
     * \ingroup building_blocks
     */

    class Population : public PopBase
    {
    public:
        //! Default constructor, creates an empty population
        Population();
        ~Population();
        //! Specify a random number generator
        void RandomNumberGenerator(RandomNumberBase *);
        //! Reconstruct population
        void Reconstruct();
        //! Enforce periodic boundary conditions on all chains
        inline void EnforceBC() {for (int i=0;i<nc;++i) chv[i]->EnforceBC();}

        int calcNSpecies();
        //@{
        /**
          @name Adding molecules to the system
          */
        //! Clear all chains
        void clear();
        //! Add proteins to the population
        /**
         * @param ntg Name of the N terminal capping group like "Acetyl"
         * @param sq The amino acid sequence, like "GEWTYDDATKTFTVTE"
         * @param ctg Name of the C terminal capping group like "Amide"
         * @param hwmny Number of chains of the specified kind you want to add.
         * It is alright to say "none" for the capping groups.
         * It is alright to add several copies of one peptides and then several
         * copies of another.
         */
        int AddProtein(std::string ntg,std::string sq, std::string ctg,
                       int hwmny=1);

        //! Add protein sequences to the population from a PDB file
        /**
         * @param hwmny Number of copies of the sequence to be added
         * @param lst A list of selected residues. The selections sould for
         * instance, come from a PDB file using mk_selection function in
         * PDBReader. SelRes objects contain chain information. So, if
         * more than one chain is detected in the list, more than one
         * chain will be added. If further, hwmny is greater than 1,
         * each chain in lst will be added hwmny times.
         */
        int AddProtein(std::list<SelRes> &lst, int hwmny=1);
        //! Add chains from a PDB file.
        /**
        * This is provided only for backward compatibility. The function
        * AddProtein(int hwmny, std::list<SelRes> &lst) above should be preferred.
        */
        int AddProtein(int hwmny, std::string pdbfilename);
        //! Add hwmany chains of a sequence described by fullseq
        /**
        * Introduced in version 1.1.0. The sequence description in fullseq
        * includes the N- and C- terminal capping groups if they should
        * be included. By default, the string is interpreted word for
        * word, each word being translated into a residue or a capping
        * group. It does not matter if you use single letter or 3-letter
        * codes or full names in those words, if the words are separated
        * by spaces. The read-mode toggles to-and-from character-mode
        * if the "*" character is encountered. In the character-mode, each
        * letter is interpreted as a one-letter symbol for an amino-acid.
        * Examples:<br>
        * fullseq="ALA ALA ALA" means alanine-alanine-alanine
        * fullseq="ALA *ALA* ALA" means alanine-alanine-leucine-alanine-alanine
        *
        * This function is useful if there is no good one-letter symbol
        * for a group, like Acetyl, D-proline etc.
        */
        int AddProtein(std::string fullseq, int hwmany=1);

        //! Assign only sequence info from an XML node
        int assign_sequences(prf_xml::XML_Node *pnode);

        //! Set "cis"-peptide-bond between residue iaa and iaa+1 in chain ich
        inline void setCis(int ich, int iaa) {
            cisres[ich].push_back(iaa);
        }

        //! Set up whether (un)charged chain ends are to be used
        /**
        * If this function is called with a false for N or C terminus, any chain
        * for which no explicit end group is specified, gets a "VoidEG" for that
        * end group. If an end group is specified, that is used. Note that
        * "uncharged chain ends" does not mean NH2 at the N-terminus and COOH
        * at the C-terminus. It just means that the terminal amino acids are
        * created just like any other, with no extra atoms.
        *
        * To use charged chain ends for un-capped sequences, this function does
        * not ever need to be called. That is the default behaviour.
        */
        inline void charged_ends(bool b1, bool b2) {charged_NTs=b1;charged_CTs=b2;}
        //@}

        //@{
        /**
          @name Assigning 3D structure to members of the population
          */
        //! Import structure from a list of AtomRecords
        /**
         * Takes a list of AtomRecords, \c rec, possibly exported by a
         * PDBReader or the Population at another time.
         *
         * This function assigns coordinates given in a list of AtomRecords
         * to the atoms in a popultion. It is useful to think of it as a
         * list copy operation. The population is like a list, and the contents
         * (coordinates) of another list (list of AtomRecords) is imported.
         * The naming of chains in the list of AtomRecords is used only
         * to separate blocks meant for different chains, i.e., the actual
         * names of the chains are ignored. The chain specified by \c at_chain
         * (default value = 0) is used as the target of the first chain in the
         * AtomRecord list.
         *
         * The argument \c assignments is a pre-allocated array of bool
         * which is used to store which atoms were actually assigned to. It
         * should be initialized elsewhere, so that it has the same size as
         * the number of atoms in the population, and all entries should be
         * initialized to false. Entries corresponding to atoms, which receive
         * new coordinates through this function, are changed to "true". The
         * other elements of the array are not touched, so that this function
         * can be called many times to assign to different parts of the
         * population. The final values in the \c assignments array can
         * be used to infer all the atoms which were assigned to.
         */
        int ImportStructure(std::list<AtomRecord> &rec,
                            std::vector<bool> &assignments, int at_chain=0);
        //! Try to infer missing coordinates
        /**
         * The argument \c assignments specifies which atoms have been
         * assigned coordinates and which atoms not. This function tries to
         * guess where those atoms with unknown coordinates should be put.
         * In reality, there is not much action in this matter in the
         * Population class. Here there is a loop over chains and the
         * corresponding function for each chain is invoked. Guessing
         * unspecified coordinates now only takes into account the known
         * geometry of protein chains, and not the non-bonded interactions.
         */
        int guess_missing_coordinates(std::vector<bool> &assignments);
        //@}

        //@{
        /**
          @name Reading in internal coordinates
          */
        //! Aggressively assign structure from an XML Node
        /**
          * If the XML node contains more chain objects than are currently
          * present in the population, new chains will be added. Then each
          * chain is forced to adopt the sequence of the corresponding chain
          * in the XML node. After this, the internal coordinates are read
          * from the XML node and assigned to the protein chains.
        */
        int Read_XML(prf_xml::XML_Node *pnode);
        //! Assign coordinates from an XML node
        /**
          * Population can be assigned a structure from an XML node, for
          * instance, as a part of the initialisation. The XML node must have
          * a name \c population, and it must have a few special child
          * tags. There could be a series of child tags of name \c protein
          * with a node structure like in ProFASi's XML output structures.
          * In addition, one can make assignments to any degree of freedom.
          * This example should be clear enough:
          * \verbatim
            <population>
              <dof_assignments>
                <dof id="::b:25"> 2.38972</dof>
                <dof id="::b:26"> -2.77682</dof>
                ...
              </dof_assignments>
            </population>
          * \endverbatim
          * The DOF id is is a string identifier for a degree of freedom. The
          syntax is described in \ref dof_strings .
          */
        int assign_structures(prf_xml::XML_Node *pnode);
        //! Read compressed binary configuration data
        inline void ReadConf(FILE *fp) {
            for (int i=0;i<nc;++i) chv[i]->ReadConf(fp);
        }
        //! Read raw configuration data in plain text format
        inline void ReadConf_text(FILE *fp) {
            for (int i=0;i<nc;++i) chv[i]->ReadConf_text(fp);
        }
        //@}

        //@{
        /**
          @name Initializing the population
          */
        //! Allocate memory and create protein objects
        /**
         * Initialize, by default creates the peptides with random values for
         * all degrees of freedom. That's how they are created, and it is
         * normally the desired starting condition in a simulation.
         *
         * From version 1.0.1, one can optionally pass an argument "1", to
         * create all proteins in the population in a "stretched out" state.
         * In case there is more than one protein chain in the system, the
         * relative position of chains will still be random. Further options
         * to start from possible cristalline geometries in multi-chain
         * systems are under consideration, and may be provided in the future
         * for different values of the optional argument.
         *
         * If a totally different starting condition is required, it can be
         * arranged after the call to Initialize. Use the Chain(i) function
         * to get a pointer to one chain. Then initialize each chain in
         * whichever way you want. Finally, if you wish randomize the relative
         * locations of different chains with the RandomizeRelConf series
         * of functions.
         */
        void Initialize(int inittyp=0);
        //! Initialize coordinates with type specified by a string
        /**
        * Introduced in version 1.1.0. The argument init_type is a
        * description of the initialization. It could be have the
        * following values: <br>
        * <ul>
        * <li> "random" : random values to all degrees of freedom </li>
        * <li> "stretched" : streched chains.</li>
        * <li> "stretched random_rel" : stretched chains, but with random relative
        * positions and orientations.</li>
        * <li> "file:somefile" : read in degrees of freedom from a text configuration
        * file </li>
        * </ul>
        */
        int InitCoord(std::string init_type);
        int Init();
        int re_index();
        int index_dof();
        //! Check consistency of the DOF index
        int check_DOF_index();

        inline bool initialized() { return initzd; }
        //! Random values to all degrees of freedom, and reconstruct system.
        void Randomize();
        //! Randomize leaving internal coordinates untouched
        void RandomizeRelConf();
        //! Randomize by moving the chain number ich rigidly
        void RandomizeRelConf(int ich);
        //! Randomize by moving the chains from ich to jch rigidly
        void RandomizeRelConf(int ich,int jch);
        //! Randomize only the internal coordinates
        void RandomizeIntConf();
        //! Randomize only the internal coordinates of chain ich
        void RandomizeIntConf(int ich);
        //! Randomize only the internal coordinates of chains ich to jch
        void RandomizeIntConf(int ich,int jch);
        //@}
        //@{
        /**
          @name Accessing constituents
          */
        //! Access i'th protein chain through a pointer
        inline Protein * Chain(int i) {return (chv[i]);}

        //! Access the longest sequence in the system by a pointer
        inline Protein * LongestChain() {return chv[longest];}

        //! Access the protein with the shortest sequence in the system
        inline Protein * ShortestChain() {return chv[shortest];}

        //! i'th ligand in the system, including all proteins, capping groups..
        inline Ligand * ligand(int i) {return lgv[i];}

        //! i'th amino acid in the system, including all protein chains
        inline AminoAcid *amino_acid(int i) {return aav[i];}

        //! Name or sequence of i'th protein chain
        inline std::string PepName(int i) {return chv[i]->Sequence();}

        //! Number of different species of Proteins in the system.
        inline int NSpecies() {return naam;}

        //! A copy of the i'th atom in the system.
        inline Atom atom(int i) const {return Atom(i,atyp[i]);}

        //! Atom type information for the i'th atom
        inline AtomKind SpeciesOf(int i) const {return atyp[i];}

        //! Total number of chains
        inline int NumberOfChains() const {return nc;}

        //! Total number of amino acids in all chains together
        inline int NumberOfResidues() const {return (int)aav.size();}

        //! Total number of ligands in all chains together
        inline int NumberOfLigands() const {return (int)lgv.size();}

        //! Total number of atoms
        inline int NumberOfAtoms() const {return ntotatms;}

        //! Number of residues in the i'th chain
        /**
        * Overrides num_res function from PopBase.
        */
        int num_grp(int i) const ;

        //! Global index of first ligand of a chain
        /**
        * Returns the integer index of the first ligand of i'th chain in
        * the vector of all ligands.
        */
        inline int chain_start(int ich) {
            return chv[ich]->memberLigand(0)->UniqueId();
        }

        //! Index of one past the last ligand of i'th chain
        inline int chain_end(int ich) {
            return 1+
                    chv[ich]->memberLigand(
                            chv[ich]->numLigands()-1)->UniqueId();
        }

        std::string chain_name(int i) const ;

        int index_of_grp(std::string ires,int ich);
        std::string grp_name(int ires, int ic);
        Ligand * existing_group(int ires, int ic);
        //@}

        //@{
        /**
          @name Managing degrees of freedom
          */
        //! Get the i'th DOF in the system
        double get_dof(size_t i);
        //! Set DOF i'th DOF value
        void set_dof(size_t i, double vl);

        //! Get info on DOF with index i in the entire system
        inline DOF_Info & get_dof_info(size_t i) {return dof_info[i];}
        //! Get info on DOF with index i within chain ich
        DOF_Info & get_dof_info(size_t ich, size_t i);

        //! Get DOF value using a DOF_Info object as key
        double get_dof(DOF_Info &d);
        //! Set DOF value using a DOF_Info object as key
        void set_dof(DOF_Info &d, double vl);

        //! Retrieve all "degrees of freedom" in a single array
        /**
          * The degrees of freedom contain all torsional DOF from all chains.
          * In addition, there is (slightly redundant) information on the
          * rigid body coordinates. 6 DOF per chain would be sufficient. But
          * reconstructing chains from such rigid body coordinates involves
          * more steps than the redundant coordinates used in PROFASI, where
          * the cartessian coordinates of the first 3 atoms of every backbone
          * are stored: 9, instead of 6 rigid body coordinates. The layout in
          * the array is coordinates of one chain followed by the other: i.e.,
          * information about one chain appears contiguously.
          */
        void get_dof(std::vector<double> &dofary);

        //! Set all "degrees of freedom" from a given array
        /**
          * See clarification on "degrees of freedom" in the documentation of
          * get_dof(std::vector<double> &dofary) above. The size of the array
          * has to be correct. No checks are performed.
          */
        void set_dof(std::vector<double> &dorary);

        //! Number of coordinates from which the exact state can be restored
        inline int n_dof() {return ndof;}

        //! Reference to the map (vector of DOF_Info) of all DOF indexes
        inline std::vector<DOF_Info> & dof_map() { return dof_info; }
        //! Interpret a string as a DOF identifier
        /**
          This function maps a ProFASi DOF identifier string to a unique
          integer global index for that degree of freedom. If the DOF can not
          be interpreted within the current population, -1 will be returned.
          \sa dof_strings
          */
        int get_dof_id(std::string dofstr);
        //! Set DOF by interpreting DOF identifier and value from strings
        void set_dof(std::string dofstr, double vl);

        //@}

        //@{
        /**
          @name Writing structure in various formats
          */
        //! Write population in XML, pdb, binary or text conf format
        /**
          * The state of a population can be written in many formats. There
          * is the PDB format. But there are other formats preserving more
          * information about the configuration. PROFASI has 3 such formats.
          * The text and binary conf formats are trivial records of the
          * degrees of freedom of one chain after the other. All numbers
          * written are "double" values. These values can be read in later
          * by the \em same population. The binary and text conf formats
          * do not contain information about what chains were present in the
          * population when the configuration was written out.
          *
          * The preferred format is XML. It is more compact than the PDB
          * format, as the PROFASI XMl files contain only the degrees of
          * freedom (torsion angles), like the binary and textconf formats.
          * But unlike those two, the XML format keeps sequence information,
          * and is a self-contained record of the population. A population
          * can be initialized from scratch using such a snapshot.
          *
          * This function handles writing in all the above mentioned formats.
          *
          * @param in_format : 1 means PDB, 2 means XML, 3 means textconf,
          *     4 means binary conf and 0 means write nothing.
          * @param flnm: Name of the snapshot file
          * @param ittime : Some "time info", typically number of MC sweeps
          * @param tindex : A "temperature index"
          * @param en: Energy.
          * Note that the Population does not know anything about temperature,
          * energy, and has no concept of any kind of MC time. It is useful
          * to have such info in the snapshots, but such info must be provided
          * to the population from outside. For backward compatibility, we do
          * not write the MC time and temperature index in the text and
          * binary configuration files created with this function.
          */
        void SaveSnapshot(int in_format, std::string flnm,
                          unsigned long ittime, int tindex, double en);
        //! Write down all the proteins in plain text
        void Write();
        //! Write some information about all the proteins
        void WriteShort();
        //! Write into binary configuration file
        inline void WriteConf(FILE *fp) {
            for (int i=0;i<nc;++i) chv[i]->WriteConf(fp);
        }

        inline void WriteConf_text(FILE *fp) {
            for (int i=0;i<nc;++i) chv[i]->WriteConf_text(fp);
        }

        //! Write population info in an XML format
        /**
        * The XML format contains both sequence and structure information for
        * the chains. The population node only contains the number of chains,
        * and a bunch of child nodes corresponding to the chains.
        */
        void Write_XML(FILE *op);

        //! Make an XML node object containing information on the population
        prf_xml::XML_Node * make_xml_node();

        //! Write PDB header lines (SEQRES and such lines before ATOM lines)
        void writePDBHeader(FILE *fp, unsigned long itime, int tindex, double entot);
        void writeSequenceInfo(FILE *fp);
        //! Export PDB to file specified by a FILE pointer
        void WritePDB(FILE *fp);
        //! Export PDB with heavy atoms first for each amino acid
        void WritePDB2(FILE *fp);

        int descriptors(std::list<SelRes> &slc,
                        std::list<AtomDescriptor> &des);
        //! Append the PDB Atom descriptor information to the end of the list
        int export_descriptors(std::list<AtomDescriptor> &lst);
        //! Make a Shape object out of the coordinates of specified atoms
        int export_shape(std::vector<int> &vct, Shape &shp);
        //@}



        inline int PairType(int i,int j) const {
            return (5*((int)atyp[i])+(int)atyp[j]);
        }
    private:
        int ntotatms,longest,shortest,naam,ndof;
        bool initzd,charged_NTs, charged_CTs;
        std::vector<Protein *> chv;
        std::vector<Ligand *> lgv;
        std::vector<AminoAcid *> aav;
        std::vector<AtomKind> atyp;
        std::vector<DOF_Info> dof_info;
        DOF_Index dofindex;
        std::deque<std::string> tmpchv;
        std::map<int, std::vector<int> > cisres;
        RandomNumberBase *rnd;
        RandomNumberBase triv;
    };
}
/**
  \page dof_strings ProFASi DOF identifier strings

  \section specs Specification
ProFASi uses the following syntax to label a degree of freedom:

\li A DOF string has 4 fields, separated by the ':' character
\li A DOF string contains no spaces
\li The first field is the chain number, an integer starting from 0
\li The second field is the residue number, an integer starting from 0.
In case there are N- or C-terminal capping groups for a chain, they
are counted as residues.
\li The third field is a type specifier for the DOF. Following values
are possible:
<ul>
   <li> "rigid_body_xyz" or "rigid" or "r" or "0"</li>
   <li> "backbone_torsional_angle" or "backbone" or "b" or "1"</li>
   <li> "sidechain_torsional_angle" or "sidechain" or "s" or "2"</li>
</ul>
\li The last field is an index which changes meaning depending on
the other fields. It represents the serial number of a DOF in the
context defined by the other fields.
\li The chain, residue and type fields may or may not be given. The
index field is compulsory.
\li If the type field is given, the index field refers to a serial
number of the DOF within its own type. "43rd backbone DOF" etc.
Otherwise, the index field refers to a serial number for the DOF
in a list of all kinds of DOF.
\li If the chain field is empty, the residue field as well as the
index field is interpreted relative to the whole system.
\li If the chain field is given, but residue field is empty, the
index is interpreted as the serial number in the given chain starting
from the beginning of the chain.
\li If both chain and residue fields are given, the residue field is
interpreted within that chain. The index is local to the residue and
within the type of the DOF.
\li It is not allowed to specify the residue field but not
specify the DOF type.

\section examples Examples
\li <b>:::22</b> is the 23rd degree of freedom in the whole population. The
chain, residue and type fields are not specified. Therefore the index 23 has the
widest possible context.
\li <b>: :b:22</b> is the 23rd backbone degree of freedom counting all chains.
The type of the DOF is specified, but the chain and residue are not given. So,
the context is "all backbone degrees of freedom in the system".
\li <b>3:::22</b> is the 23rd degree of freedom in the 4'th chain (indices
start from 0).
\li <b>3: :s:7</b> is the 8th side chain degree of freedom in the 4'th chain
\li <b>:41:s:2</b> is the 3rd side chain DOF of the 42nd residue in the whole
system. N- and C-terminal capping groups count as residues. The residue index
is relative to the whole system.
\li <b>0:12:s:2</b> is the 3rd side chain DOF of the 13th residue of the first
chain. Notice that whenever the residue field is specified, the type field is
also specified. This is a requirement.
  */
#endif

