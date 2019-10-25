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

#ifndef MINIMIZER_HH
#define MINIMIZER_HH
#include "../Aux/NDFunctional.hh"
#include "../Aux/LineProjection.hh"

class Minimizer
{
public:
    Minimizer();
    ~Minimizer();
    void set_function(NDFunctional *gf);
    inline void set_reset_freq(int n) {ncgreset=n;}
    inline void set_max_cyc(int n) {maxcyc=n;}
    inline void set_tolerance(double v) {ftol=v;}
    inline void set_update_strategy(int i) {updttyp=(i==0)?0:1;}
    double minimize();
    inline void get_pos(std::valarray<double> &v) {f->get_position(v);}
    inline void set_pos(std::valarray<double> &v) {f->set_position(v);}
    inline void set_good_enough(double x) {
        goodenough=x;
        knowgoodenough=true;
    }
    inline void set_abs_scale(bool vl=true) {absscale=vl;}
    inline void unset_good_enough()
    {
        knowgoodenough=false;
    }
private:
    NDFunctional *f;
    LineProjection lp;
    int ndof,ncgreset,maxcyc,updttyp;
    double ftol,goodenough;
    bool knowgoodenough,absscale;
};

#endif // MINIMIZER_HH
