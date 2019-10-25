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

#ifndef ProteinRMSD_HH
#define ProteinRMSD_HH
#include "../Elements/Protein.hh"
#include "../Aux/rmsd.hh"
#include "../Aux/PDBReader.hh"
#include "../Aux/RMSD_Utils.hh"
#include "Observable.hh"

namespace prf
{
    //! RMSD evaluator between a file and a collection of atoms in the program
    /**
     * \ingroup profasi_observables
     * This is an "Observable" class which facilitates RMSD evaluation from a
     * structure in a PDB file with a structure in the program. Neither the
     * structrue in the file nor in the program need to be complete Proteins.
     *
     * The two comparison structures are built in a flexible way. One can
     * select non-overlapping regions in the PDB file as follows:<br><br>
     * ProteinRMSD myrmsd;<br>
     * myrmsd.set_struc1("abc.pdb:7:A,12,19:A,25,45");<br><br>
     * "abc.pdb" is the name of the reference pdb file. The MODEL 7 in that file
     * is selected. In that model, the residues 12-19 and 25-45 of chain A are
     * selected. The set of all selected atoms constitutes the comparison structure.
     * Any number of segments can be selected, and they need not be from the
     * same chain. But they have to be from the same MODEL. The model id can be
     * omitted, so that it defaults to the first model. If residue range is not
     * specified, the entire chain is selected. <br>
     *
     * The convention for residue range specification here is different from
     * the rest of PROFASI. The range 12-19 above means the residues labeled
     * as 12 through 19 (inclusive at both ends) in the file. In case the residue
     * numbering in the file jumps from 13 to 16, the selection would include
     * residues 12,13,16,17,18,19. It is precisely to deal with pdb files with
     * omitted residues or with insertion codes, that the residue number labels
     * are merely treated as strings.
     *
     * A corresponding comparison object can be created from the PROFASI Population.
     * The syntax is exactly the same, but instead of the file name, one uses the
     * symbol $. For example, <br>
     * myrmsd.set_struc2("$::A,12,19:A,25,45");<br><br>
     * Notice that the numbering here starts from 1, to keep it consistent with
     * the selection numbering in pdb files.
     *
     * It is possible to set both struc1 and struc2 to be simulated segments in the
     * program, instead of some external pdb file. For instance, one might be
     * interested in the similarity between the corresponding segments in
     * chain A and chain C in the program, without reference to any external pdb
     * file, like...
     *
     * myrmsd.set_struc1("$::A,1,15");<br><br>
     * myrmsd.set_struc2("$::C,1,15");<br><br>
     *
     * One does not always want to include all the atoms in the selected segments
     * for an RMSD calculations. To filter only the interesting ones, there are,
     * well, filters! To use the backbone atoms and CB atoms only, you would write,
     *
     * myrmsd.filter("+BB+CB"); <br>
     * If you additionally want to exclude proline residues and sulfur atoms from
     * the calculations, then you would write instead,<br>
     *
     * myrmsd.filter("+BB+CB-%PRO-@S");<br>
     *
     * It is assumed by default that the first of the two structures does
     * not change during the simulation. This is normally the case, like
     * when it is an NMR derived PDB file.  Some calculations related to
     * that structure are done only once. If you use the population as the
     * source of the first structure, by default, it would mean the state
     * of the population at the time of invocation of ProteinRMSD::init.
     * If on the other hand, you wish a genuine changing state to be compared
     * with another changing state, you have to declare the first state to
     * be mobile, using the following function.
     *
     * myrmsd.mobile_first_struc();
     *
     *
     * Finally, you need the RMSD value for whatever selections and filters you
     * made. This is done in the usual way, as in all other classes derived from
     * Observable, call refresh(), and then Value(); //available through inheritance.
     *
     * myrmsd.refresh();
     * x=myrmsd.Value();
     */

    class ProteinRMSD : public Observable
    {
    public:
        ProteinRMSD(); //!< Default constructor
        virtual ~ProteinRMSD();
        int init_obs();
        //! Set the first structure, see the discussion at the class description
        inline void set_struc1(std::string rf) {fl1=rf;}

        //! Set the second structure
        inline void set_struc2(std::string rf) {fl2=rf;}

        //! Set default selection
        inline void set_default_sel1(std::string s1) {dsel1=s1;}
        inline void set_default_sel2(std::string s1) {dsel2=s1;}

        //! Declare the first structure as changable
        inline void live_first_struc(bool tf=true) { fixed1=!tf; }

        //! Add a filter for the kind of atoms to be used in the RMSD calculation
        void filter(std::string flt);
        //! Don't try to align sequences
        /**
        * With this option, the program does not try to align sequences. It is
        * useful, in combination with appropriate filters, to calculate RMSD
        * between similar but not identical sequences.
        */
        void no_seq_align() {seq_align=false;}

        //! If called with "true", the observable value becomes exp(-MSD/100)
        /**
         * The quantity Q=exp(-RMSD^2 /100) has a range between 0 and 1 and
         * has been found useful. We include an option for the Value() function
         * to return this.
         */
        inline void returnQ(bool tf) {calcQ=tf;}

        //! Estimate range for RMSD
        virtual void rangeEstimate(double &x1, double &x2);
        //! Initialize with selected ranges and filter
        int init();
        //! Delete information for the second object
        /**
        * This should be called if you wish to call init over and over
        * again, with new files for the second structure.
        */
        inline void delete_matrix2() {pdb2.delete_matrix();}

        double evaluate();
        double gradientXYZ(std::valarray<double> &g);
        void set_logger_threshold(int thr);
    protected:
        void make_mapped_lists(std::list<AtomDescriptor> &l1,
                               std::list<AtomDescriptor> &l2,
                               std::list<AtomDescriptor> &ml1,
                               std::list<AtomDescriptor> &ml2);
        PopBase *obj1,*obj2;
        std::string fl1,fl2,sel1,sel2,myfilter,altcrd,dsel1,dsel2;
        std::vector<int> sel_atoms1,sel_atoms2;
        PDBReader pdb1,pdb2;
        Shape shape1,shape2;
        RMSD rmsd;
        bool fixed1,calcQ, seq_align;
        double rnge;
        RMSD_Utils utils;
    };
}

#endif
