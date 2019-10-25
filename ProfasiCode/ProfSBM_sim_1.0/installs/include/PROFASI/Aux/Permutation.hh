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

#ifndef Permutation_HH
#define Permutation_HH
#include "profasi_io.hh"
#include <vector>

namespace prf_utils
{
    //! Permutation class originally developed by SM for ALFS
    /**
     * Permutation class taken from ALFS: Area Law Fragmentation of Strings.
     * ALFS is a particle physics Monte Carlo event generator that handles
     * hadronization of a given partonic state, based on the Lund model area
     * law. In connection with symmetrisation of amplitudes for final states
     * with many identical mesons, several utility classes were developed,
     * of which Permutation and ExchangeLink are found to be convenient in
     * connection with parallel tempering. ExchangeLink has been renamed to
     * IntPair, because we think it is a better more desciptive name. These
     * classes are otherwise only slightly modified to conform to the profasi
     * conventions, and prf_utils namespace.
     *
     */

    class Permutation
    {
    public:
        Permutation();
        //! Create a permutation of n objects
        Permutation(int n);
        ~Permutation();
        //! assignment
        Permutation & operator=(Permutation & prm);
        //! Whether this is an identity permutation
        bool IsIdentity();
        //! checks that index i is in the appropriate range
        bool Element(int i);
        //! Content of the i'th slot, no range check
        inline int operator[](int i) {return valperm[i];}

        //! Slot location for the integer i
        inline int Transpose(int i) {return trperm[i];}

        //! Change size to i
        void reset(int i);
        //! Restore to the identity permutation
        void SetToIdentity();
        //! Assign val to the location indx
        int Assign(int indx, int val);
        //! Exchange objects at positions i and j
        /**
         * Exchange objects at positions i and j, irrespective of what they
         * are. This means that if the current state of the permutation of 3
         * integers is {b,a,c}, and you call FlipContents(1,3) you get
         * {c,a,b}, i.e., contents of positions i and j are exchanged.
         */
        int FlipContents(int i, int j);
        //! Exchange the integers i and j, wherever they are
        /**
         * Exchange objects a and b, irrespective of where they
         * are. This means that if the current state of the permutation of 3
         * integers is {b,a,c}, and you call FlipLocations(a,b)
         * you get {a,b,c}, i.e., the locations of the objects represented by
         * integers a and b, are exchanged, so that a ends up where b was,
         * and vice versa.
         */
        int FlipLocations(int a, int b);
        friend prf::Output & operator<<(prf::Output & os, Permutation &prm);
    private:
        int nmembers;
        std::vector<int> valperm;
        std::vector<int> trperm;
    };
}

#endif
