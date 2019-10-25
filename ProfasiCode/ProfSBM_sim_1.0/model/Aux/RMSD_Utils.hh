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

#ifndef RMSD_Utils_HH
#define RMSD_Utils_HH
#include "CompositePredicate.hh"
#include <list>
#include <vector>

namespace prf
{

    class AtomDescriptor;

    class SelRes;

    //! A trivial predicate to select all atoms
    /**
      * \ingroup pdb_handling
      */

    class isAtom : public SimplePredicate<AtomDescriptor>
    {
    public:
        isAtom();
        ~isAtom();
        bool operator()(AtomDescriptor atd);
    };

    //! A predicate that returns true if the atom is not Hydrogen
    /**
    * \ingroup pdb_handling
    */

    class isHeavyAtom : public SimplePredicate<AtomDescriptor>
    {
    public:
        isHeavyAtom();
        ~isHeavyAtom();
        bool operator()(AtomDescriptor atd);
    };

    //! A predicate that returns true if the atom is on the backbone
    /**
    * \ingroup pdb_handling
    */

    class isBackboneAtom : public SimplePredicate<AtomDescriptor>
    {
    public:
        isBackboneAtom();
        ~isBackboneAtom();
        bool operator()(AtomDescriptor atd);
    };

    //! A predicate that returns true if the atom is of a certain species
    /**
    * Useful if you want to select all carbon atoms, or exclude only
    * sulfur atoms in a selection.
    * \ingroup pdb_handling
    */

    class AtomTypeSelector : public SimplePredicate<AtomDescriptor>
    {
    public:
        explicit AtomTypeSelector(char gtyp);
        ~AtomTypeSelector();
        bool operator()(AtomDescriptor atd);
        char mytp;
    };

    //! A predicate that returns true if the atom is in one kind of residue
    /**
    * One can invert it to use it as a filter against, for instance,
    * Proline residues in an RMSD calculation, if needed.
    * \ingroup pdb_handling
    */

    class ResSelector : public SimplePredicate<AtomDescriptor>
    {
    public:
        explicit ResSelector(std::string gtyp);
        ~ResSelector();
        bool operator()(AtomDescriptor atd);
        std::string mytp;
    };

    //! A predicate that returns true if the atom has a certain label
    /**
    * Could be set up to return true if the label of an atom, as found in its
    * AtomDescriptor object, matches the filter value. By default, if set to
    * select atoms of label " CG ", it will also accept CG1 and CG2 as matches.
    * But it can be told to be strict, and only return true upon an exact match.
    * \ingroup pdb_handling
    */

    class LabelSelector : public SimplePredicate<AtomDescriptor>
    {
    public:
        explicit LabelSelector(std::string lbl);
        ~LabelSelector();
        bool operator()(AtomDescriptor atd);
        //! Use strict label matching
        void set_strict(bool x) {strct=x;}

        std::string label;
        bool strct;
    };

    //! The RMSD_Utils handler class.
    /**
    * A helper class to perform tasks useful in connection, for instance,
    * with the evaluation of RMSD. It converts a string in the PROFASI filter
    * format (see \ref filters) to an actual filter. It can take two lists of
    * AtomDescriptors and create a mapping between the atoms, deleting unmatched
    * entries from the given lists.
    * \ingroup pdb_handling
    */

    class RMSD_Utils
    {
    public:
        typedef CompositePredicate<AtomDescriptor> CPA;
        typedef SimplePredicate<AtomDescriptor> SPA;
        RMSD_Utils();
        ~RMSD_Utils();
        int make_filter(std::string expr);
        void add_filter(std::string flt, std::list<CPA > &);
        inline std::string filter_name() {return total_filter.id();}

        void clear_filters();
        void apply_filter(std::list<AtomDescriptor> &lst);
        int map_atoms(std::list<AtomDescriptor> l1,std::list<AtomDescriptor> l2,
                      std::vector<int> &v1,std::vector<int> &v2 ,
                      std::string opts="");
        void alt_coord_policy(std::string pol);
        void alt_coord_parse(std::list<AtomDescriptor> &lst);
        bool seq_match(std::list<SelRes> &l1,std::list<SelRes> &l2);
        void seq_align(std::list<SelRes> &l1,std::list<SelRes> &l2);
        int score(std::vector<std::string> &v1, std::vector<std::string> &v2, size_t i1, size_t i2);
        void align_segs(std::list<SelRes>::iterator i1, std::list<SelRes>::iterator i2,
                                std::list<SelRes>::iterator j1, std::list<SelRes>::iterator j2,
                                std::list<SelRes> &resi, std::list<SelRes> &resj);
        bool clean_seq(std::list<SelRes> &reslist);
        inline void set_logger_threshold(int vl) {log_thres=vl;}
    private:
        std::list<CPA> incl_list,excl_list,cumul_p,cumul_m;
        std::list<SPA *> predlst;
        CPA total_filter;
        std::string mypolicy;
        int log_thres;
    };
}

#endif
