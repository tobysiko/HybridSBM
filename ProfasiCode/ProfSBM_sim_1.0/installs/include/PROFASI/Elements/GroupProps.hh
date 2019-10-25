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

#ifndef GroupProps_HH
#define GroupProps_HH
#include <string>
#include <vector>
#include <map>
#include "AtomKind.hh"

namespace prf
{

    class GroupProps
    {
    public:
        GroupProps();
        GroupProps(std::string lbls);
        ~GroupProps();
        void init(std::string lbls);
        void set_type(std::string typ);
        void torsion_dof(int ai[4]);
        void torsion_dof(int a1, int a2, int a3, int a4);
        inline bool isAA() const {return isaa;}

        inline bool isEG() const {return iseg;}

        inline bool isSynth() const {return issynth;}

        inline bool isNat() const {return isnat;}

        inline bool isHydrophobic() const {return hyph;}

        inline bool isPolar() const {return !hyph;}

        inline bool isCharged() const {return chgd;}

        inline int num_atoms() const {return natoms;}

        inline int num_hvatoms() const {return nhvatoms;}

        inline int num_sdof() const {return nsdof;}

        inline std::string label(int i) const {return thelabel[i];}

        inline int index(std::string lbl) {return theindex[lbl]-1;}

        bool valid_label(std::string lbl) const;

        inline std::vector<std::string> valid_labels() const { return thelabel; }

        AtomKind species(int i);
        void dof_def_atoms(int idof, int &a1,int &a2,int &a3,int &a4);
        void dof_def_atoms(int idof, int ai[4]);
        inline std::string CommonName() const {return comname;}

        inline void CommonName(std::string st) {comname=st;}

        inline std::string TLC() const {return threelet;}

        inline void TLC(std::string st) {threelet=st;}

    private:
        bool isaa,iseg,issynth,isnat,hyph,chgd;
        int natoms,nhvatoms,nsdof;
        std::vector<std::string> thelabel;
        std::map<std::string,int> theindex;
        std::vector<int> defatoms;
        std::string comname,threelet,typedescr;
    };
}

#endif
