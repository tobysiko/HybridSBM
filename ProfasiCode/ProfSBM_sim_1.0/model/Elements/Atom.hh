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

#ifndef Atom_HH
#define Atom_HH
#include <vector>
#include "AtomCoordinates.hh"
#include "../Aux/profasi_io.hh"
#include "AtomKind.hh"

/**
* \defgroup building_blocks Building Blocks
* The building blocks are a collection of classes representing specific
* concepts we use when we think of a protein system. These classes try
* to translate these concepts into code, which is then used by other
* "method" classes to perform a calculation.
*/

/**
  * \defgroup physical_entities Physical entities
  * \ingroup building_blocks
  * @brief Representation of physical entities relevant for ProFASi
  *
  * This module groups a set of ProFASi classes representing things like
  * the atoms, the amino acids, proteins etc.
  */

namespace prf
{
    //! Representation of an atom in PROFASI
    /**
     * The Atom is a very simple class with one member defining which kind of
     * atom it is, and another specifying its coordinates in space. The
     * coordinates are represented by an AtomCoordinates object. For many
     * purposes, it is sufficient to think of the position as an ordinary
     * three vector. An atom has a so called "unique_id" which can be
     * accessed through the function UniqueId().
     * \ingroup physical_entities
     */

    class Atom
    {
    private:
        AtomKind my_type;
        AtomCoordinates pos;
    public:
        static std::vector<double> radius;
        Atom();
        Atom(int i);
        Atom(int i, AtomKind tp);
        Atom(const Atom &atm);
        ~Atom();
        Atom & operator=(const Atom &ga);   //!< Clones both type and position
        //! The unique integer id for the atom
        inline int UniqueId() const {return pos.Id();}

        inline void UniqueId(int i) {pos.Id(i);}

        //inline void Register() {pos.atyp[pos.Id()]=my_type;}
        //! Access position as a three vector
        inline Vector3 Position() const {return pos.value();}

        //! Access position as an AtomCoordinates object.
        inline AtomCoordinates Pos() const {return pos;}

        //! Assign a position to the atom
        inline void Pos(const Vector3 &gpos) {pos=gpos;}

        inline void translate(const Vector3 &trv) {pos+=trv;}

        //! Return the type of atom, H, C, N...
        inline AtomKind Species() const {return my_type;}

        inline void Species(AtomKind i) {my_type=i;}

        //! Returns its name as a C string
        const char *strName() const;
        void Write();
        friend Output & operator<<(Output &op, const Atom &a);
    };
}

#endif
