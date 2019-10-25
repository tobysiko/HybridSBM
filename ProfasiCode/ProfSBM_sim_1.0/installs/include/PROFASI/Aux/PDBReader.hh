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

#ifndef PDBReader_HH
#define PDBReader_HH
#include <string>
#include <vector>
#include <list>
#include "../Elements/PopBase.hh"
#include "profasi_io.hh"
#include "AtomRecord.hh"
#include "Shape.hh"
#include "Matrix.hh"
#include "fileutils.hh"

namespace prf
{
    //! An abstraction of one model in a PDB file
    /**
    * When a PDB file is parsed, it is convenient to find out what models
    * it contains, and store enough information about the models, so that
    * a particular model can be easily selected for further processing.
    * \ingroup pdb_handling
    */

    class PDB_model
    {
    public:
        PDB_model();
        ~PDB_model();
        PDB_model(const PDB_model &);
        PDB_model & operator=(const PDB_model &);
        //! Translate residue number written in file
        /**
        * Residue numbers written in a PDB file may start from 634, may
        * not even be strictly numeric in nature. This function takes one
        * such number label, and translates it into an integer index,
        * measured from the first residue in the file for that chain.
        *
        * The first residue is numbered 0, no matter what it is called in
        * the PDB file. Given any PDB file, however bizarre, there exists
        * this numbering relative to the first residue in that file, which
        * is well defined. At least that is the faith.
        */
        int index_of_res(std::string chid, std::string resnum);
        //! Return index of a given chain id
        int index_of(std::string chid);
        int num_chains;
        std::vector<std::string> chain_name;
        std::vector<std::vector<std::string> > res_map;
        std::vector<std::vector<int> > res_ln0,res_lnn;
        void print_summary(prf::Output &op);
    };

    //! A top level PDB file parser
    /**
    * The PDBReader processes a PDB file, selects a (not necessarily contiguous)
    * region of residues in the file, and can return either a list of
    * AtomRecords or AtomDescriptors corresponding to the selection.
    * \ingroup pdb_handling
    */

    class PDBReader : public PopBase
    {
    public:
        PDBReader();
        //! Constructor with an implicit call to set_file
        PDBReader(std::string fl);
        ~PDBReader();

        //! Sets up PDBReader to work with a given file
        void set_file(std::string fl);
        //! Select model i for further processing
        void set_model(int i);
        //! Return currently used file name
        inline std::string file_name() const {return myfile;}

        //! Read in the PDB file as a character matrix
        /**
        * Expensive operation. Returns 2 and quits without doing anything if it has
        * already been called on a file. A new file can be read only after calling
        * delete_matrix(). Returns 1 if it succeeds reading the file into a character
        * matrix. Returns 0 if it fails for some reason. Allocates the matrix as
        * required. Implicitly calls gather_model_info()
        * \sa delete_matrix(), gather_model_info()
        */
        int read_matrix();
        //! Discards stored character matrix of the PDB file.
        /**
        * In the interest of good memory management, when it is clear that the program
        * wont need to make a fresh selection on a PDB file, it is good to release
        * the large amount of memory taken by the character matrix.
        */
        void delete_matrix();
        int clear();

        //! Atom descriptor information for a given SelRes list
        /**
        * For a given list of selected residues, put  the
        * AtomDescriptor objects corresponding to all atoms in those
        * residues in a given list. Contents of the list before this
        * function is called will be lost.
        */
        int descriptors(std::list<SelRes> &slc,
                        std::list<AtomDescriptor> &lst);
        //! AtomRecord objects for the selection list
        /**
        * One can export AtomRecord objects as well, which includes
        * both the AtomDescriptors and coordinates.
        */
        int records(std::list<SelRes> &slc, std::list<AtomRecord> &lst);

        //! Make a Shape object out of the coordinates of specified atoms
        /**
        * One of the most important functions that really makes this
        * class usable. Once selections have been exported, and the
        * exported list has been subject to various filtres, one might
        * want to know the 3-dimensional shape of the filtred selection.
        * In a way, this enables very complicated selections for, say,
        * RMSD evaluations. One can first select a range in the PDB file.
        * Then export the atom descriptors. Then apply filtres so that only
        * a desired subset of atoms remain. Then use this function to get
        * a Shape object consisting only of these selected atoms. The
        * integers in the vector<int> argument must specify the indices
        * internal to this class, which can be used to retrieve a certain
        * AtomDescriptor or AtomRecord. When exporting the descriptors or
        * records, such an index is imprinted in the atom descriptor objects.
        * So, when external filtres have been applied, one can retrieve
        * the necessary integer indices in the remaining AtomDescriptors.
        * The relevant public member of the AtomDescriptor class is
        * called int_label. The vector vct passed to this function is expected
        * contain those internal identifiers for the atoms as given by
        * this class.
        *
        */
        int export_shape(std::vector<int> &vct, Shape &shp);
        //! Similar to export_shape, but returns the records instead
        int export_records(std::vector<int> &vct,
                           std::list<AtomRecord> &shp);
        //! Similar to export_shape, but returns only AtomDescriptors
        int export_descriptors(std::vector<int> &vct,
                               std::list<AtomDescriptor> &shp);
        std::string res_index(int ar_indx,std::string chid);
        //!Sequence of chain ich as a list of 3 letter codes
        int sequence_of_chain(int ich,std::list<std::string> &seq);
        //! Label of the i'th chain in the current model
        std::string chain_name(int i) const;
        //! Number of groups (residues and capping groups) in chain ic
        int num_grp(int ic) const;
        //! Natural index of a group with a string index ires in chain ich
        /**
        * The "residue serial number" in PDB files is, in general, best not
        * thought of as a number. Sometimes it will contain non-digit
        * characters. Sometimes, even if all residue serial numbers are integers
        * the difference between two of them may not represent the actual number
        * of residues between them. This is because in PDB files there are
        * often missing residues. The "natural index" is the index relative to
        * what is in the file. It is always a natural number.
        *
        * As an example, if a PDB file starts with the first few residues numbered
        * 41,42,43,44, the natural index of the one numbered 44 is 3. The counting
        * starts from 0.
        */
        int index_of_grp(std::string ires, int ich);
        //! Residue or capping group name for residue with natural index ires
        std::string grp_name(int ires,int ich);
        //! The string index corresponding to the residue with natural index ires
        std::string str_index(int ires,int ich);
        //! Return reference to the AtomLabelDictionary used to interpret labels
        inline AtomLabelDictionary & dictionary() {return aldict;}
        bool check_file();
    private:
        void gather_model_info();
        std::string read_word(int rwo,int colm,int wdth);
        void write_word(int rwo, int colm,std::string txt);
        void mk_rec(int lno, AtomRecord &r);
        Vector3 mk_vec(int lno);
        int add_records(SelRes &sel, std::list<AtomRecord> &rcd);
        int add_descriptors(SelRes &sel, std::list<AtomDescriptor> &des);
        bool check_unknown_labels();
        std::string myfile;
        std::vector<PDB_model> models;
        int nto;
        Matrix<char> data;
        bool matrix_read;
        AtomLabelDictionary aldict;
    };
}

#endif
