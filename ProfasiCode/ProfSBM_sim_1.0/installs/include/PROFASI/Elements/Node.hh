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

#ifndef Node_HH
#define Node_HH
#include "Atom.hh"
#include <vector>
#include <deque>
#include "../Aux/ConnectionsMatrix.hh"
#include "../Aux/Named.hh"

/**
* \defgroup geometrical_objects Geometrical objects
* The geometrical objects coded in ProFASi are those which recur frequently
* in protein structures. Bonds connect atoms, and atoms can be thought of as
* junctions of bonds. These junctions could have 4 or 3 legs, the legs could
* be symmetrically or asymmetrically placed and so on. But the total number of
* essentially different geometrical constructs needed to put together a protein
* chain is relatively small.
* \ingroup building_blocks
*/

namespace prf
{
    //! A Node is (like in graph theory) a meeting point of a few bonds
    /**
     * The Node class is used to facilitate implementation of different kinds
     * molecules. A Node has one incoming bond, and a few out going bonds.
     * There is of course one atom sitting at the node junction, and one atom
     * attached to the incoming and each of the out going bonds. All these
     * atoms are regarded as members of a node. By definition, a Node can have
     * at most one degree of freedom: rotation about its incoming bond. This
     * degree of freedom can sometimes be frozen or locked, if the molecule
     * can not rotate about the incoming bond of the junction. The
     * outgoing bonds can themselves lead to other nodes, which will be called
     * subnodes of the parent node.
     *
     * The atom located at the junction of bonds will be referred to as
     * "Junction". The atom connected to it by the incoming bond is the "Root".
     * Atoms connected to the Root (other than Junction) can serve as the
     * reference direction to define a torsional angle about the incoming
     * bond. Such an atom is called a "Base". There is one primary Base and
     * upto two additional "alternative" Bases.
     * \ingroup geometrical_objects
     */

    class Node : public Named
    {
    public:
        Node();
        Node(const Node &);
        Node & operator=(const Node &);
        virtual ~Node();
        //@{
        /**
         * @name Atom assignment to nodes
         * The AssignAtoms series of functions is a means of allowing an
         * an external routine to associate a node with any pre-defined set of
         * atoms. Exactly how the atoms are to be organized inside a node
         * to facilitate calculations is left for the inheritted classes to
         * implement, according to the nature of the involved node. By default,
         * it is assumed that you pass the Base,Root,Junction followed by a
         * reference to an external atom array with an index from which to start
         * reading the rest of the atoms. This behaviour is sufficient for
         * simple short node types where the atoms in the node are simply
         * covalently connected to the junction atom.
         */

        virtual void AssignAtoms(Atom &a0, Atom &a1,Atom &a2,
                                 std::vector<Atom> &att,int st);
        virtual void AssignAtoms(Atom &a0, Atom &a1,Atom &a2,Atom &a3);
        virtual void AssignAtoms(Atom &a0, Atom &a1,Atom &a2,
                                 Atom &a3,Atom &a4);
        virtual void AssignAtoms(Atom &a0, Atom &a1,Atom &a2,
                                 Atom &a3,Atom &a4,Atom &a5);
        //@}
        //! Initialize the node, possibly pre-calculate anything that can be
        virtual void Initialize();
        //! Assign coordinates to atoms attached to outgoing bonds
        /**
         *Using the position of Base, Root and Junction, and the state of the
         * degree of freedom associated with the node, find the coordinates
         * of all other atoms in the node.
         */
        virtual void Create();
        void ReCreate(); //!< Recursive Create. Self followed by subnodes
        virtual void AssignPhi(double);   //!< assign to the dof of the Node.
        inline double Phi() {return phi;} //!< retrieve the torsional angle

        //! Lock the torsional dof of the Node to value xx
        inline void LockPhi(double xx) {phi=xx;rotlocked=1;}

        inline int Locked() const {return rotlocked;}//!< if it is really a dof

        inline Atom & ATOM(int i) {return atm[i];}   //!< Reference to i'th atom

        inline Atom & atom(int i) {return atm[i];}   //!< Reference to i'th atom

        void AddSubnode(Node *nd);   //!< Make another node a subnode
        void AtomOffset(int i);   //!< Offset UniqueIds of all atoms by i
        //! Revert the torsional angle to a given value.
        /**
         * Unlike AssignPhi, RevertPhi does not implicitly calculate the
         * coordinates of all attached atoms. This function is meant to be
         * called when an Update is rejected, and assumes that the atom
         * coordinates have been restored from the backup values externally.
         */
        void RevertPhi(double);
        //@{
        /**
         * @name Bond length assignments
         * The following four functions refer to the maximum number of such
         * parameters that one might need, assuming there are at most 3
         * out-going bonds.
         */
        virtual void SetBondLengths(double, double, double, double, double);
        virtual void SetBondAngles(double,double,double,double);
        virtual void SetRelPhi(double,double);
        virtual void SetBranchLength(int i, double x);
        //@}
        //! Export connection information to an external ConnectionsMatrix
        virtual void ExportConnections(ConnectionsMatrix &aa);
        //! List pairs of member atoms which are connected by 3 covalent bonds
        virtual void LocPairs(std::deque<std::pair<int,int> > & lcp);
        //! Assign proper Bases to subnodes and call BuildConnections on them
        virtual void BuildConnections();
        //! Set a new atom as one Base
        virtual void SetBase(Atom & ap);
        //!be careful: This actually means set Base, Root and Junction
        void SetBase(Atom &,Atom &,Atom &);
        //!And this means set all bases
        /**
        * Any previously stored bases will be forgotten, and only the values
        * given here will be used.
        */
        void set_bases(Atom &b1, Atom &b2, Atom &b3);
        void set_bases(Atom &b1, Atom &b2);
        //! Atoms that move when the node turns
        /**
         * Returns an inclusive range between r1 and r2. This means all atoms
         * with unique_ids between r1 and r2 (inclusive) should be assumed to
         * have moved.
         */
        void MobileAtoms(int & r1, int &r2);
        //! Set what range should be returned as the mobile atoms.
        void SetMobileAtoms(int i, int j);
    protected:
        double phi;
        int rotlocked;
        std::vector<Atom> atm;
        Atom altbase[3];
        int naltbases;
        Node *subnode[2];
        int nsubnd;
        int upst,upnd;
    };
}

#endif
