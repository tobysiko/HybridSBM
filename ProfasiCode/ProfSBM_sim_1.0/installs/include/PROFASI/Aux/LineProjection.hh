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

#ifndef LINEPROJECTION_HH
#define LINEPROJECTION_HH
#include "UnaryFunction.hh"
#include "NDFunctional.hh"

class LineProjection : public UnaryFunction
{
public:
    LineProjection();
    ~LineProjection();
    void set_base_function(NDFunctional *gf);
    void set_origin(std::valarray<double> &org);
    void set_direction(std::valarray<double> &drc);
    inline void get_last_position(std::valarray<double> &vl) {vl=fpos;}
    double operator()(double t);
protected:
    NDFunctional *ndf;
    std::valarray<double> orig,dirc,fpos;
};

#endif // LINEPROJECTION_HH
