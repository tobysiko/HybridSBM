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

#ifndef RandomPermutation_HH
#define RandomPermutation_HH
#include "RandomNumberBase.hh"
#include <vector>

namespace prf_utils
{
    //! Create a random permutation of objects in a given vector
    /**
    * This function uses random numbers given in a separate vector (ran)
    * to permute the contents of a given vector of arbitrary type (res).
    */
    template <class T>
    void random_permutation(std::vector<T> &res, std::vector<double> &ran)
    {
        std::vector<T> img=res;
        int remains=res.size(),indx=0;

        for (size_t i=0;i<res.size();++i) {
            indx=(int) remains*ran[i];
            res[i]=img[indx];
            img[indx]=img[--remains];
        }
    }

    //! Create a random permutation of objects in a given vector
    /**
    * This function uses random numbers generated by a given random number generator
    * to permute the contents of a given vector of arbitrary type.
    */
    template <class T>
    void random_permutation(std::vector<T> &res, prf::RandomNumberBase *ran)
    {
        std::vector<double> v(res.size(),0);

        for (size_t i=0;i<res.size();++i) v[i]=ran->shoot();

        random_permutation<T>(res,v);
    }
}

#endif
