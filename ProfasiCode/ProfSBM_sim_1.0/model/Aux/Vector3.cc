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

#ifndef VECTOR3_CC
#define VECTOR3_CC
#include "Vector3.hh"

namespace prf
{
    Vector3::Vector3(double tx, double ty, double tz)
            : dx(tx), dy(ty), dz(tz) {}

    Vector3::Vector3(const Vector3 & p)
            : dx(p.x()), dy(p.y()), dz(p.z()) {}

    Vector3::~Vector3() {}


    void Vector3::theta(double th)
    {
        double ma   = mag();
        double ph   = phi();

        x(ma*sin(th)*cos(ph));
        y(ma*sin(th)*sin(ph));
        z(ma*cos(th));
    }

    void Vector3::phi(double ph)
    {
        double ma   = mag();
        double th   = theta();

        x(ma*sin(th)*cos(ph));
        y(ma*sin(th)*sin(ph));
        z(ma*cos(th));
    }

    void Vector3::mag(double ma)
    {
        double th = theta();
        double ph = phi();

        x(ma*sin(th)*cos(ph));
        y(ma*sin(th)*sin(ph));
        z(ma*cos(th));
    }

    double Vector3::angle(const Vector3 & q) const
    {
        double ptot2 = mag2()*q.mag2();
        return ptot2 <= 0.0 ? 0.0 : std::acos(dot(q)/sqrt(ptot2));
    }

    void Vector3::rotateX(double gangle)
    {
        double s = sin(gangle);
        double c = cos(gangle);
        double yy = y();
        y(c*yy - s*z());
        z(s*yy + c*z());
    }

    void Vector3::rotateY(double gangle)
    {
        double s = sin(gangle);
        double c = cos(gangle);
        double zz = z();
        z(c*zz - s*x());
        x(s*zz + c*x());
    }

    void Vector3::rotateZ(double gangle)
    {
        double s = sin(gangle);
        double c = cos(gangle);
        double xx = x();
        x(c*xx - s*y());
        y(s*xx + c*y());
    }

    void Vector3::rotate(double gangle, const Vector3 & axis)
    {
        rotateZ(-axis.phi());
        rotateY(-axis.theta());
        rotateZ(gangle);
        rotateY(axis.theta());
        rotateZ(axis.phi());
    }

    Vector3 operator * (double a, const Vector3 & p)
    {
        return Vector3(a*p.x(), a*p.y(), a*p.z());
    }

    double Vector3::torsion(const Vector3 &prev, const Vector3 &next)
    {
        Vector3 v12,v23,v212;
        v12=-cross(prev);
        v23=cross(next);
        v212=cross(v12);
        return std::atan2(v12.mag()*v23.dot(v212),v212.mag()*v23.dot(v12));
    }
}

#endif

