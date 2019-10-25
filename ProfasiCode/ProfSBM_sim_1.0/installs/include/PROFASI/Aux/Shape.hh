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

#ifndef Shape_HH
#define Shape_HH
#include <vector>
#include "profasi_io.hh"
#include "Vector3.hh"
#include "Matrix.hh"

namespace prf
{
    //! Represents an abstract shape in space defined by a set of points
    /**
     * The Shape is a helper class for the similarity measure RMSD. A shape
     * in space is defined by the coordinates of an enumerated set of points.
     * A Shape can be tranlated and rotated in space, and in general acted
     * upon by a linear transformation.
     */

    class Shape
    {
    public:
        typedef std::vector<Vector3> PointRepresentation;
        Shape(); //!< Create an empty Shape object
        virtual ~Shape();
        Shape(unsigned int i);   //!< Shape object with i points
        void Resize(unsigned int i);   //!< Change number of points to i
        //! Specify the location of i'th point
        inline void Point(int i, const Vector3 &v) {pnts[i]=v;}

        //! Get the location of the i'th point
        inline Vector3 Point(int i) const {return pnts[i];}

        //! Return the whole vector of points
        PointRepresentation & Points() {return pnts;}

        inline int NPoints() const {return nref;} //!< Number of Points

        void Translate(const Vector3 &t);   //!< translate shape by vector t
        //! Rotate shape by angle th about axis ax, through center of mass
        void Rotate(double th,const Vector3 &ax);
        //! Scale x, y and z directions with different factors
        void Scale(double xscl, double yscl, double zscl);
        //! Linear transformation by applying a matrix R on every vector
        /**
         * Assumes that for mtR, the coordinates are column vectors
         */
        void Apply1(Matrix<double> & mtR);
        //! More linear transformation by applying a matrix R on every vector
        /**
         * Assumes that for mtR, the coordinates are row vectors
         */
        void Apply2(Matrix<double> & mtR);
        //! Save shape in a file
        virtual void Save(const char *fln);
        //! find center of mass
        void find_cm();
        //! Find CMS and translate coordinates such that it is at origin
        void find_and_move_to_cm();
        //! Get center of mass
        inline Vector3 center_of_mass() const {return cm;}

        //! Get radius of gyration
        inline double gyr_r2() const {return meansq;}

    protected:
        int nref;
        Vector3 cm;
        double meansq;
        PointRepresentation pnts;
    };
}

#endif
