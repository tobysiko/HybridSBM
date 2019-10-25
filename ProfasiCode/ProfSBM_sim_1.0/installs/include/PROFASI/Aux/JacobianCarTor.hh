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

#ifndef JACOBIANCARTOR_HH
#define JACOBIANCARTOR_HH
#include "Matrix.hh"
#include "../Elements/Population.hh"

// Jacobian matrix of Cartessian coordinates with respect to the torsion angles
class JacobianCarTor : public Matrix<double>
{
public:
    JacobianCarTor();
    ~JacobianCarTor();
    void set_population(prf::Population * popu);
    void set_dofmap(std::vector<int> &indxs);
    void refresh();
private:
    prf::Population *popl;
    std::vector<std::pair<int, int> > axis, chlist;
    std::vector<int> chstrt,indexmap;
    int nat;
};

#endif // JACOBIANCARTOR_HH
