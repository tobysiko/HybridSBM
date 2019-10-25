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

#ifndef DOF_INDEX_HH
#define DOF_INDEX_HH
#include <vector>
#include <map>
#include <vector>
#include "DOF_Info.hh"

namespace prf {
    class key_pair
    {
    public:
        key_pair();
        key_pair(int,int);
        ~key_pair();
        key_pair(const key_pair &);
        key_pair &operator=(const key_pair &);
        bool operator==(const key_pair &) const;
        bool operator<(const key_pair &) const;
        inline void assign(int i1, int i2) { k1=i1;k2=i2; }
        inline int key1() const { return k1; }
        inline int key2() const { return k2; }
    private:
        int k1,k2;
    };
    class key_tripplet
    {
    public:
        key_tripplet();
        key_tripplet(int,int,int);
        ~key_tripplet();
        key_tripplet(const key_tripplet &);
        key_tripplet &operator=(const key_tripplet &);
        bool operator==(const key_tripplet &) const;
        bool operator<(const key_tripplet &) const;
        inline void assign(int i1,int i2,int i3) {k1=i1;k2=i2;k3=i3;}
        inline int key1() const { return k1; }
        inline int key2() const { return k2; }
        inline int key3() const { return k3; }
    private:
        int k1,k2,k3;
    };
    //! Helper class for population for DOF management
    /**
      This class provides additional structure to the vector of degrees of
      freedom in the population, so that the unique id of a DOF can be
      querried about a partially known DOF_Info object. The DOF_Info class
      stores redundant information on degrees of freedom, so that a few
      combinations of its fields suffice to identify it uniquely. For instance,
      specifying the chain, the DOF type and the index of the DOF in its type
      is sufficient to uniquely identify it. The same DOF might also be found
      by specifying the chain, the residue, the type and the type specific
      index in the residue. When initialized the DOF_Index class creates a
      map of all complete descriptions of a DOF on to its unique id.
      */
    class DOF_Index
    {
    public:
        DOF_Index();
        ~DOF_Index();
        //! Reset all previous info
        void reset();
        //! Initialize with respect to a DOF array
        int init(std::vector<DOF_Info> &dofs);
        //! Get the unique id from chain and index in chain
        int get_uid_from_chain_index(int ch, int chind);
        //! Get the unique id from type and the global index within type
        int get_uid_from_type_index(int tp, int tpind);
        //! Get the unique id from chain, type and chain level index in type
        int get_uid_from_chain_type(int ch, int tp, int chtpind);
        //! Get the unique id from residue, type and residue level index in type
        int get_uid_from_residue_type(int rs, int tp, int rstpind);

    private:
        std::map<key_pair, int> ch_chind, tp_tpind;
        std::map<key_tripplet, int> ch_t_chtind, rs_t_rstind;
    };
}
#endif // DOF_INDEX_HH
