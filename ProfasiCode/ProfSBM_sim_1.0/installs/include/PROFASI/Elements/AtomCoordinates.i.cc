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

using namespace prf;

AtomCoordinates & AtomCoordinates::operator= (const Vector3 &g)
{
    (*uc)[locn]=g.x();
    (*uc)[locn+1]=g.y();
    (*uc)[locn+2]=g.z();
    return *this;
}

Vector3 AtomCoordinates::operator+ (const AtomCoordinates &g) const
{
    return Vector3((*uc)[locn]+(*uc)[g.locn],
                   (*uc)[locn+1]+(*uc)[g.locn+1],
                   (*uc)[locn+2]+(*uc)[g.locn+2]);
}

Vector3 AtomCoordinates::operator+ (const Vector3 &g) const
{
    return Vector3((*uc)[locn]+g.x(),
                   (*uc)[locn+1]+g.y(),
                   (*uc)[locn+2]+g.z());
}

AtomCoordinates & AtomCoordinates::operator+= (const AtomCoordinates &g)
{
    int gloc=g.locn;
    (*uc)[locn]+=(*uc)[gloc];
    (*uc)[locn+1]+=(*uc)[gloc+1];
    (*uc)[locn+2]+=(*uc)[gloc+2];
    return *this;
}

AtomCoordinates & AtomCoordinates::operator+= (const Vector3 &g)
{
    (*uc)[locn]+=g.x();
    (*uc)[locn+1]+=g.y();
    (*uc)[locn+2]+=g.z();
    return *this;
}

Vector3 AtomCoordinates::operator- (const AtomCoordinates &g) const
{
    int gloc=g.locn;
    return Vector3((*uc)[locn]-(*uc)[gloc],
                   (*uc)[locn+1]-(*uc)[gloc+1],
                   (*uc)[locn+2]-(*uc)[gloc+2]);
}

Vector3 AtomCoordinates::operator- (const Vector3 &g) const
{
    return Vector3((*uc)[locn]-g.x(),
                   (*uc)[locn+1]-g.y(),
                   (*uc)[locn+2]-g.z());
}

AtomCoordinates & AtomCoordinates::operator-= (const AtomCoordinates &g)
{
    int gloc=g.locn;
    (*uc)[locn]-=(*uc)[gloc];
    (*uc)[locn+1]-=(*uc)[gloc+1];
    (*uc)[locn+2]-=(*uc)[gloc+2];
    return *this;
}

AtomCoordinates & AtomCoordinates::operator-= (const Vector3 &g)
{
    (*uc)[locn]-=g.x();
    (*uc)[locn+1]-=g.y();
    (*uc)[locn+2]-=g.z();
    return *this;
}

double AtomCoordinates::dot(const AtomCoordinates &g) const
{
    int gloc=g.locn;
    return (*uc)[locn]*(*uc)[gloc]+
           (*uc)[locn+1]*(*uc)[gloc+1]+
           (*uc)[locn+2]*(*uc)[gloc+2];
}

double AtomCoordinates::dot(const Vector3 &g) const
{
    return value().dot(g);
}

Vector3 AtomCoordinates::operator*(const AtomCoordinates &g) const
{
    return value()*g.value();
}

Vector3 AtomCoordinates::operator*(const Vector3 &g) const
{
    return value()*g;
}

double AtomCoordinates::dist2(const int i1, const int i2)
{
    return sqr(xdist(i1,i2)) +sqr(ydist(i1,i2)) +sqr(zdist(i1,i2));
}

double AtomCoordinates::seps2(const int i, const double xpt, const double ypt,
                              const double zpt)
{
    return sqr(sepx(i,xpt)) +sqr(sepy(i,ypt)) +sqr(sepz(i,zpt));
}

double AtomCoordinates::pseps2(const int i, const double xpt, const double ypt,
                               const double zpt)
{
    return sqr(psepx(i,xpt)) +
           sqr(psepy(i,ypt)) +
           sqr(psepz(i,zpt));
}
