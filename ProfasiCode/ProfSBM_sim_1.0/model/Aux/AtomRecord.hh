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

#ifndef AtomRecord_HH
#define AtomRecord_HH
#include <string>
#include "Vector3.hh"
#include "AtomLabelDictionary.hh"

/**
* \defgroup pdb_handling ProFASi PDB handling module
* This is a collection of classes that implements reading and
* parsing of PDB files in ProFASi. The top level class, PDBReader
* is probably the only one most users are likely to encounter
* directly. It reads a PDB file. It can be told to select one
* area of the file, and apply certain filters on it, and then
* "export" the resulting list of AtomRecords, or AtomDescriptors.
*
* At a lower level, this module contains a scheme for making
* arbitrary combination of filters, which have been used in this
* module to include or exclude atoms based on any property an
* atom record in a PDB file might have. But the filter composition
* methods here have potentially a wider applicability.
\ingroup utilities
\sa \ref mimiqa
*/

namespace prf
{
    //! Constants describing the format of an atom record in PDB files
    /**
    * \ingroup pdb_handling
    */

    namespace prf_pdb_vars {
        const int chain_label_column = 21; //!< Chain label start column
        const int chain_label_width = 1; //!< Chain label width
        const int res_num_column=23; //!< Residue number start column
        const int res_num_width=3; //!< Residue number width
        const int res_label_column=17; //!< Residue label start column
        const int res_label_width=3; //!< Residue label width
        const int atom_label_column=12; //!< Atom label column
        const int atom_label_width=4; //!< Atom label width
        const int atom_index_column=6; //!< Atom index column
        const int atom_index_width=5; //!< Atom index width
        const int atom_type_column=13; //!< Atom type column
        const int alt_coord_column=16; //!< Alternative records column
        const int crd_x_col=30; //!< X coordinate start column
        const int crd_y_col=38; //!< Y coordinate start column
        const int crd_z_col=46; //!< Z coordinate start column
        const int crd_col_width=8; //!< Width of coordinates fields
        const int occup_col=56; //!< Occupancy column
        const int occup_width=4; //!< Width of occupancy field
    }

    //! Minimal information that identifies one atom in a PDB file
    /**
    * To identify one particular atom line in a PDB file, one needs
    * a residue number, an atom label. Certain other properties are
    * fixed by these two for a given molecule: like the residue name
    * species of the atom etc. The AtomDescriptor class is an
    * abstraction for such an incomplete description of the atoms.
    * The description is only incomplete in the sense that it does not
    * specify any coordinates for the atom. This is useful in many
    * circumstances where the coordinates are irrelevant.
    * \ingroup pdb_handling
    */

    class AtomDescriptor
    {
    public:

        //! nickname
        /**
        * int_label is a label given to the descriptor from an outside program or
        * function. It is not a property of the atom represented by the
        * descriptor. It is like a "nickname".
        */
        int int_label;
        //! model number
        int imdl;
        //! atom index
        int  iatom;
        //! Relative residue index
        /**
        iresrel is the residue serial number starting from the first residue
        * present in the file, counting only the residues explicitly present
        * in the file.
        */
        int iresrel;
        //! Alternative coordinates identifier
        /**
        * ialt is a character label to distinguish between alternative coordinates
        * for one atom, sometimes given in a PDB file. Most often this character
        * will be blank. When there are two alternative coordinates given,
        * they will differ in their values of ialt.
        */
        char ialt;
        //! Atom species
        /**
        * atom_type is H,C,N, O or S depending on if the atom is hydrogen, carbon,
        * nitrogen, oxygen or sulfur, respectively
        */
        char atom_type;
        //! Identifies the type of the record represented by the PDB line.
        /**
        * keyword is the string formed by the first 6 characters in a PDB line.
        * It designates what kind of record that line is. For an atom, its value
        * is either "ATOM  " or "HETATM". resnm is the residue name in 3 letter
        * code.
        */
        std::string keyword;
        //! 3 letter code for the residue
        std::string resnm;
        //! Atom identifier string inside a residue
        /**
        * atom_label is the label identifying the atom inside one residue,
        * like " CA ", " OG2"  etc.
        */
        std::string atom_label;
        //! ich is the chain label.
        std::string ich;
        //! ires is the residue index as given in the PDB file.
        /**
        * It is treated as a string and not as a number, because that is how it is
        * often treated in PDB files. A sequence of consecutive residue numbers in
        * a PDB file could read 81, 81A, 81B, 82 ... Or, it could be 15, 16, 24, 25, ...
        * There is no clean way to treat this index as an integer and cover all cases.
        * One thing we assume is that  if the residues are different, the index will
        * at least be different!
        */
        std::string ires;
        //! occ is the occupancy value often mentioned for Atoms in PDB files.
        std::string occ;
    public:
        AtomDescriptor();
        ~AtomDescriptor();
        AtomDescriptor(const AtomDescriptor &);
        //! Assignment operator
        AtomDescriptor &operator=(const AtomDescriptor &);
        //! Comparison of two AtomDescriptors
        /**
        * This is a correspondence function rather than an equality test. The
        * optns argument could take values "ignore_res_index", "ignore_res_label",
        * "ignore_chain_label", or "ignore_model_label" to make the comparison
        * insensitive to a particular property. We may want the comparison to be
        * insensitive to the residue label when we want to find backbone RMSD
        * between two proteins of different sequence, for instance.
        */
        bool corresponds_to(AtomDescriptor &,std::string optns="");
        //! Set fields from a string, assumed to be a PDB line
        int set_fields(std::string, AtomLabelDictionary *d);
        //! Imprint fields onto a 80 character wide line, intended to be a PDB line
        void mark_fields(char gline[81]);   //81, not 80. There is the \0 as well.
        //! Try to make a good ProFASi style atom label out of a given string
        /**
         * This function maps some commonly found labels for atoms, in
         * particular hydrogen atoms, to the ProFASi convention for their
         * labels. For instance, " HA2" and " HA3" are mapped to "1HA " and " 2HA"
         * respectively.
         */
        void filter(std::string &alabel);
        //! Print
        void print_info();
        //! Return a brief description as a string, for use in output elsewhere.
        std::string short_info();
    };

    //! Record corresponding to one atom in a PDB file
    /**
    * Essentially, this is an AtomDescriptor plus 3D coordinates for an atom.
    * \ingroup pdb_handling
    */

    class AtomRecord
    {
    private:
        //! All information in a PDB atom record except coordinates
        AtomDescriptor dsc;
        //! 3D coordinates of the atoms.
        Vector3 coord;
        std::string line;
    public:
        AtomRecord();
        AtomRecord(std::string, AtomLabelDictionary *d);
        ~AtomRecord();
        //! Set up fields based on a PDB line
        int set_fields(std::string, AtomLabelDictionary *d);
        //! Create a PDB line from stored information
        void build_pdb_line_from_fields();
        //! Access to the contained AtomDescriptor object
        inline AtomDescriptor & descriptor() { return dsc; }

        //! Access to the coordinates
        inline Vector3 coordinates() { return coord; }

        //! The corresponding PDB line
        inline std::string pdb_line() { return line; }

        //! Comparison with another record
        inline bool operator==(AtomRecord &rcd) {return line==rcd.line;}

        //! Print
        void print_info();
    };
}

#endif
