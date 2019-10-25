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

#include "Dipole.hh"

namespace prf
{
    Dipole::Dipole() :i1(0),i2(0), st(1), chid(0), ligid(0) {}

    Dipole::~Dipole() {}

    Dipole::Dipole(int gi1,int gi2) :i1(gi1),i2(gi2),
            st(1), chid(0), ligid(0)  {}

    Dipole::Dipole(int gi1,int gi2,double gst) : i1(gi1),i2(gi2),st(gst)
    {
        ligid=chid=0;
    }

    Dipole::Dipole(const Dipole &dp)
    {
        i1=dp.i1;
        i2=dp.i2;
        st=dp.st;
        chid=dp.chid;
        ligid=dp.ligid;
    }

    Dipole & Dipole::operator=(const Dipole &dp)
    {
        if (this!=&dp) {
            i1=dp.i1;
            i2=dp.i2;
            st=dp.st;
            chid=dp.chid;
            ligid=dp.ligid;
        }

        return *this;
    }

    void Dipole::SetAtoms(int gi1,int gi2) {i1=gi1;i2=gi2;}

    void Dipole::SetStrength(double gst) {st=gst;}

}
