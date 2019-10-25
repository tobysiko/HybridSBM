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

#ifndef Vector3_H
#define Vector3_H
#include <cmath>

namespace prf
{
    //! A representation of three vectors
    /**
     * The code for this class originally comes from the HEP_THREE_VECTOR
     * class from the CLHEP package.  It's a small but very useful class.
     *
     * \ingroup building_blocks
    */

    class Vector3
    {
    public:
        Vector3(double x = 0.0, double y = 0.0,double z = 0.0);
        Vector3(const Vector3 &);
        ~Vector3();
        inline void set_val(double gx,double gy, double gz)
        {dx=gx;dy=gy;dz=gz;}

        //! Access and change the components in cartesian coordinates.
        //@{
        /**
         * @name Components in cartesian coordinate system.
         */
        inline double x() const {return dx;}

        inline double y() const {return dy;}

        inline double z() const {return dz;}

        inline void x(double gx) {dx=gx;}

        inline void y(double gy) {dy=gy;}

        inline void z(double gz) {dz=gz;}

        inline double operator[](int i) {
            switch (i) {
                case 0:
                    return dx;
                case 1:
                    return dy;
                default:
                    return dz;
            };
        }

        //@}
        //! The azimuth angle.
        inline double phi() const {
            return x() == 0.0 && y() == 0.0 ? 0.0 : std::atan2(y(),x());
        }

        //! The polar angle.
        inline double theta() const {
            return x() == 0.0 && y() == 0.0 && z() == 0.0 ?
                   0.0 : std::atan2(perp(),z());
        }

        //! The magnitude squared.
        inline double mag2() const {
            return x()*x() + y()*y() + z()*z();
        }

        //! The magnitude.
        inline double mag() const {return sqrt(mag2());}

        //! The transverse component squared..
        inline double perp2() const {
            return x()*x() + y()*y();
        }

        //! The transverse component.
        inline double perp() const {return sqrt(perp2());}

        //! Set phi keeping mag and theta constant.
        void phi(double);
        //! Set theta keeping mag and phi constant.
        void theta(double);
        //! Set magnitude keeping theta and phi constant.
        void mag(double);

        //! Assignment.
        inline Vector3 & operator = (const Vector3 & p) {
            x(p.x());
            y(p.y());
            z(p.z());
            return *this;
        }

        //! Addition.
        inline Vector3 & operator += (const Vector3 &p) {
            dx += p.x();
            dy += p.y();
            dz += p.z();
            return *this;
        }

        //! Addition.
        inline Vector3  operator + (const Vector3 &p) const {
            return Vector3(dx+p.dx,dy+p.dy,dz+p.dz);
        }

        //! Subtraction.
        inline Vector3 & operator -= (const Vector3 &p) {
            dx-=p.x();
            dy-=p.y();
            dz-=p.z();
            return *this;
        }

        //! Subtraction.
        inline Vector3  operator - (const Vector3 &p) const {
            return Vector3(dx-p.dx,dy-p.dy,dz-p.dz);
        }

        //! Unary minus.
        inline Vector3 operator - () const {return Vector3(-dx, -dy, -dz);}

        //! Scaling with real numbers.
        inline Vector3 & operator *= (double a) {
            dx *= a;
            dy *= a;
            dz *= a;
            return *this;
        }

        //! Scaling with real numbers.
        inline Vector3 operator * (double a) const {
            Vector3 ret = *this;
            ret *= a;
            return ret;
        }

        //! Scalar product.
        inline double dot(const Vector3 &p) const {
            return x()*p.x() + y()*p.y() + z()*p.z();
        }

        //! Cross product.
        inline Vector3 cross(const Vector3 &p) const {
            return Vector3(y()*p.z()-p.y()*z(), z()*p.x()-p.z()*x(),
                           x()*p.y()-p.x()*y());
        }

        //! Cross product.
        inline Vector3 operator*(const Vector3 &p) const {
            return Vector3(y()*p.z()-p.y()*z(), z()*p.x()-p.z()*x(),
                           x()*p.y()-p.x()*y());
        }

        //! The angle w.r.t. another 3-vector.
        double angle(const Vector3 &) const;

        //! Torsion angle about one vector
        /**
        * The vectors prev, this and next can be thought of as leading
        * from arbitrary points a->b, b->c, and c->d respectively. The
        * torsion angle is defined as the angle on the plane perpendicular
        * to the vector b->c (this), between the projections of the points.
        * The sense of the angle is positive if the rotation vector from
        * the projection of 'a' to the projection of 'd' is in the same
        * direction as the vector "this".
        */
        double torsion(const Vector3 &prev, const Vector3 &next);

        void rotateX(double);
        //!< Rotates the Vector3 around the x-axis.

        void rotateY(double);
        //!< Rotates the Vector3 around the y-axis.

        void rotateZ(double);
        //!< Rotates the Vector3 around the z-axis.

        void rotate(double, const Vector3 &);
        //!< Rotates around the axis specified by another Vector3.

        //! Scaling of 3-vectors with a real number
        friend Vector3 operator * (double a, const Vector3 &);

    private:

        double dx, dy, dz;
        // The components.
    };
}

#endif /* HEP_THREEVECTOR_H */
