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

#ifndef FindCoord_HH
#define FindCoord_HH
#include "Vector3.hh"
#include <cmath>

/**
* \defgroup utilities Utilities
* This is a diverse and loosely bound group of classes that are necessary
* for performing some actions on the concept classes grouped under the module
* \ref building_blocks .
*/

namespace prf
{
    //! Simple but very useful class used widely in PROFASI to get coordinates
    /**
     * The task performed by this class is well-defined. We have vectors u and
     * v given. We want to find a vector w, which has a certain length,
     * a certain angle with v, thvw, a certain torsional orientation about v
     * with respect to a direction defined by u.
     *
     * The vector w can be expressed as a linear combination of u, v and u*v
     * where u*v is the cross product. Often in an application, u and v would
     * have fixed lengths, angles, and even the length of w and the v-w angle
     * would be known to be fixed. So, except for the torsional angle, most
     * other constants in the linear combinations can be precalculated. This
     * is what this class does. It precalculates those constants, and for a
     * given phi angle, and vectors u, v, it quickly generates a vector w
     * satisfying all conditions.
     *
     * The phi angle here is defined with respect to the forward going
     * direction of u, so that if w were parallel to u, it would have phi=0
     *
     * Care has to be taken in using this class. The fast evaluation does
     * not check if the passed u and v vectors have the correct lengths or
     * angle between them. If they don't, the generated w would not have the
     * right properties. This means repeated, cyclic applications will lead
     * to divergence because numerical precision errors in the lengths of u
     * and w will eventually get amplified to order 1. The slow fresh_eval()
     * function will generate a vector w of required length and angles with
     * respect to v and u, no matter if the previously known lengths of u and v
     * don't match up. But it is much slower.
     *
     */

    class FindCoord
    {
    public:
        FindCoord(); //!< default constructor
        //! Lengths of u, v and w, u-v angle and v-w angle given
        FindCoord(double gu,double gv,double gw,double gthuv,double gthvw);
        //! Lengths of u, v, w, u-v and v-w angle and u-v-w torsion angle given
        FindCoord(double gu,double gv,double gw,double gthuv,double gthvw,
                  double gh);
        //! Assign properties later, if default constructor was used
        void Initialize(double gu,double gv,double gw,double gthuv,
                        double gthvw);
        //! Assign properties later, if default constructor was used
        void Initialize(double gu,double gv,double gw,double gthuv,
                        double gthvw,double gh);

        ~FindCoord();
        //! Find a vector w, given vectors u and v
        Vector3 operator()(const Vector3 &,
                           const Vector3 &);
        //! Find a vector w, given vectors u and v, and torsional angle phi
        Vector3 operator()(const Vector3 &,
                           const Vector3 &,double);
        //! Fresh evaluation of w, without using precalculated constants
        Vector3 fresh_eval(Vector3,Vector3,double);
        //! Even lazier: Assume U, V and their cross product UV are given
        Vector3 operator()(const Vector3 &,
                           const Vector3 &,const Vector3 &);
        //! U, V, U cross V and phi are given
        Vector3 operator()(const Vector3 &,
                           const Vector3 &,const Vector3 &,
                           double);
        //! Precalculate constants assuming that phi is fixed
        void FixPhi(double gph);
    private:
        FindCoord(const FindCoord &);
        double wsthvw,wcthvw;
        double sph,cph,wu,wv,wv1,wv2,wvxu;
    };
}

#endif
