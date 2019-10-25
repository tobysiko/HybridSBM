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

#ifndef NDFUNCTIONAL_HH
#define NDFUNCTIONAL_HH
#include <valarray>

class NDFunctional
{
public:
    NDFunctional();
    virtual ~NDFunctional();
    inline int n_dim() const { return nd; }
    virtual void allocate_coords(int ndim);
    virtual void set_position(std::valarray<double> &x);
    virtual void get_position(std::valarray<double> &x);
    virtual void set_coordinate(int i, double x);
    virtual double get_coordinate(int i);
    virtual double value();
    virtual double operator()(std::valarray<double> &x);
    virtual void gradient(std::valarray<double> &g);
    inline void set_scale_small(double x) {h_small=x;}
protected:
    virtual double partial_derivative(int idir);
    virtual double delta(int idir, double h);
    int nd;
    double h_small;
    std::valarray<double> cur;
};

#endif // NDFUNCTIONAL_HH
