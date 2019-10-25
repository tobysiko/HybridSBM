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

#ifndef Dipole_HH
#define Dipole_HH

namespace prf
{

    class Dipole
    {
    private:
        int i1,i2;
        //dipole moment goes from the atom placed at i1 to that at i2
        double st;
        int chid,ligid;
    public:
        Dipole();
        ~Dipole();
        Dipole(const Dipole &);
        Dipole(int,int);
        Dipole(int,int,double);
        Dipole &operator=(const Dipole &dp);
        void SetAtoms(int,int);
        void SetStrength(double);
        inline int atom1() const {return i1;}

        inline int atom2() const {return i2;}

        inline double strength() const {return st;}

        inline void Reduce(double xx) {st*=xx;}

        inline int ligand() const {return ligid;}

        inline int chain() const {return chid;}

        inline void ligand(int lg) {ligid=lg;}

        inline void chain(int ch) {chid=ch;}
    };
}

#endif
