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

#ifndef RMSD_HH
#define RMSD_HH
#include "Shape.hh"
#include "NDFunctional.hh"

namespace prf
{
    class rmsd_help_fn : public NDFunctional
    {
    public:
        rmsd_help_fn();
        ~rmsd_help_fn();
        double value();
    };

    //! Root Mean Square Deviation
    /**
     * RMSD compares two Shape objects. One Shape is set as the reference
     * structure. The other is the comparison Shape, which is passed as
     * an argument to the value functions. The root mean square deviation
     * is then defined as the minimum with respect to all possible rotations
     * and translations, of the summed squares of displacements of
     * corresponging points.
     *
     * RMSD is evaluated using closed form analytic expressions based on
     * Singular Value Decomposition.
     * \ingroup utilities
     */

    class RMSD
    {
    public:
        RMSD(); //!< Define an RMSD evaluator
        ~RMSD();
        inline double value(Shape &sh1, Shape &sh2, bool sh1incms=false) {
            return operator()(sh1,sh2,sh1incms);
        }

        //! Calculation of RMSD value between two shape objects
        /**
        * Calculates the RMSD between two shape objects. The boolean
        * variable sh1incms tells if the first shape is already in
        * the centre of mass coordinates. If not, it is translated to
        * the centre of mass coordinates during evaluation. If one of
        * the structures is fixed and is used repeatedly for rmsd
        * calculations, repeated translations of that shape to its
        * centre of mass is unnecessary, and can be avoided by
        * passing sh1incms as true. Only the first structure can be
        * fixed.
        */
        double operator()(Shape &sh1, Shape &sh2, bool sh1incms=false);
        double eval_w_grad(Shape &s1,Shape &s2,std::valarray<double> & dr);
    private:
        rmsd_help_fn f;
    };

}

#endif
