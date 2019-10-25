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

#ifndef PopBase_HH
#define PopBase_HH
#include <list>
#include "../Aux/AtomRecord.hh"
#include "../Aux/Shape.hh"

namespace prf
{
    //! An abstraction for a residue selected for further processing
    /**
    * Information, such as, which model, which chain, the residue label,
    * the residue index (possibly in a PDB file) and the true index relative
    * to the existing part of a chain representation.
    * \ingroup pdb_handling
    */

    class SelRes
    {
    public:
        SelRes();
        ~SelRes();
        SelRes(int mdlno,int chno,int i);
        SelRes(const SelRes &);
        SelRes & operator=(const SelRes &);
        bool operator==(SelRes);
        void setval(int,int,int,std::string,std::string);
    public:
        int mdl,chn,indx_nat;
        std::string indx_str,nam;
    };

    //! Base class for a population of representations of proteins
    /**
    * Two kinds of representations of proteins and their populations occur
    * inside PROFASI. One is the regular Protein and Population classes, with
    * real Atom objects whose energy can be evaluated, and which can be
    * affected by conformational Updates. But when parsing PDB files, we deal
    * frequently with the same kinds of relations between similar concepts.
    * This file is meant to represent the common features. It facilitates
    * in comparing structures between an "internal" protein object with one
    * in a PDB file. The comparison program only needs to deal with the common
    * properties.
    */

    class PopBase
    {
    protected:
        PopBase();
    public:
        virtual ~PopBase();
        //! Number of chains
        inline int num_chains() const {return nc;}

        //! Currently selected model
        inline int get_model() const {return curmodel+1;}

        //! Select model i (makes sense only for PDB files)
        virtual void set_model(int i);
        //! Label of the i'th chain
        virtual std::string chain_name(int i) const;
        //! Integer index (starting from 0) of the chain labeled chnm
        virtual int chain_number(std::string chnm) const;
        //! Number of groups (residues and capping groups)
        virtual int num_grp(int ich) const;
        //! Natural index of group labeled "ires" in the chain ich
        virtual int index_of_grp(std::string ires,int ich);
        //! String index of the chain with natural index ires
        virtual std::string str_index(int ires,int ich);
        //! Label of the group with natural index ires.
        virtual std::string grp_name(int ires,int ich);
        //! Return a list of all atom descriptors for a given selection
        virtual int descriptors(std::list<SelRes> &slc,
                                std::list<AtomDescriptor> &des);
        //! Export a pure set of 3D coordinates for particular atoms indicated by integer labels
        virtual int export_shape(std::vector<int> &vct, Shape &shp);
        int select(std::string chainno, std::string res1,
                   std::string res2,std::list<SelRes> &slc);
        //! Create a selected list of residues using a selection string slcstr
        /**
        * The string is translated into a selection criteria using \ref prf_sel_fils
        */
        int mk_selection(std::string slcstr, std::list<SelRes> &slc);
        int add_selection(int ich,int ires,std::list<SelRes> &slc);
    protected:
        int nc,curmodel,nmodels;
    };
}

#endif
