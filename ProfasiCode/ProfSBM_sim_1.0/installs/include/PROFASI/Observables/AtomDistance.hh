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

#ifndef AtomDistance_HH
#define AtomDistance_HH
#include "Observable.hh"

namespace prf
{
    //! Distance between two designated atoms
    /**
     * A somewhat more general observable than the end-to-end distance. This
     * could represent the end-to-end distance, if you just initialize it
     * with the two end atoms of the relevant backbone. But end-to-end
     * distance could just be more noisy in some situations, when there are
     * loose unstructured loops at the ends of the chain. One might then
     * want to use the distance between the structured end-points of the
     * chain. This class allows one to use any two atoms for this distance.
     * It knows nothing about atoms however, as it blindly takes the integer
     * id of the atoms, and asks AtomCoordinates class to return the distance
     * between those to ids. It is your responsibility to make sure that the
     * two indices you pass make sense.
     * \ingroup profasi_observables
     */

    class AtomDistance : public Observable
    {
    public:
        AtomDistance();
        ~AtomDistance();
        int init_obs();
        //! attach to two atoms given by their unique ids.
        inline void set_atoms(int i1,int i2) {a1=i1;a2=i2;}

        double evaluate();

        void rangeEstimate(double &xmin,double &xmax);

    private:
        int a1,a2;
    };
}

#endif
