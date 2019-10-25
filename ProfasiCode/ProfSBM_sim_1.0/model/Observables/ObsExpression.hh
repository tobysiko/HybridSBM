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

#ifndef OBS_EXPRESSION_HH
#define OBS_EXPRESSION_HH
#include "../Aux/NDFunctional.hh"
#include "../Observables/Observable.hh"
#include "../Aux/JacobianCarTor.hh"
#include "../Updates/IncrementalUpdate.hh"

class ObsExpression : public NDFunctional
{
public:
    ObsExpression();
    ~ObsExpression();
    void init();
    inline void set_population(prf::Population *popl) {
        p=popl;update.connect(popl);
    }
    void reset_dof_list();
    void reset_obs_list();
    void add_obs(double scl, Observable *o);
    inline size_t n_obs() const { return obs.size(); }
    Observable * used_obs(size_t i) {return i<obs.size()?obs[i]:NULL;}
    void add_dof(DOF_Info &d);
    double value();
    void gradient(std::valarray<double> &g);
    void get_position(std::valarray<double> &x);
    void set_position(std::valarray<double> &x);
    double get_coordinate(int i);
    void set_coordinate(int i, double x);
protected:
    double delta(int idir, double h);
private:
    Population *p;
    std::deque<Observable *> obs;
    std::deque<double> c;
    std::vector<int> indexmap,dirct,jcobi,fdelt,basic;
    std::valarray<double> gtmp,gxsum,gxtmp;
    IncrementalUpdate update;
    JacobianCarTor J;
};

#endif // ObsExpression_HH
