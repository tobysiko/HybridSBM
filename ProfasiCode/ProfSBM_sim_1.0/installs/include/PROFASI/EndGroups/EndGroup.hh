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

#ifndef EndGroup_HH
#define EndGroup_HH
#include "../AA/AminoAcid.hh"

namespace prf
{

    class EndGroup:public Ligand
    {

    protected:
        Vector3 * bv0, *bv1;
        Atom rf0, rf1;
        AminoAcid *res;
    public:
        explicit EndGroup(OneLetterCode cod);
        EndGroup(const EndGroup &);
        EndGroup & operator=(const EndGroup &);
        virtual ~ EndGroup();
        virtual bool pep_bond_link() const;
        virtual void pep_bond_atoms(int &i1,int &i2);
        virtual void Allocate();
        inline void SetRefv(Vector3 * v1, Vector3 * v2) {
            bv0 = v1;
            bv1 = v2;
        }

        inline void SetRefa(Atom & a1, Atom & a2) {
            rf0 = a1;
            rf1 = a2;
        }

        inline void Connect_to_AA(AminoAcid * acd) {
            res = acd;
        }

        virtual void Initialize();
        virtual void Reconstruct();
        virtual void Write();
        virtual void WritePDB(int &istatm, char chid, int aaind,
                              FILE * fp);
        virtual void WritePDB2(int &istatm, char chid, int aaind,
                               FILE * fp);
        virtual int rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1);
        virtual int rotDof_assign(int il, double mgd);
        virtual double get_rotDof(int il);

        virtual void LocPairsatRTdof(int i,
                                     std::deque<std::pair<int,int> >&lcp);
        virtual void BuildConnections();
        void nodes_reconnect(Atom &a1, Atom &a2);
        Node *node_for_dof(int i);
    };
}

#endif
