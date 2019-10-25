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

#ifndef INCREMENTALUPDATE_HH
#define INCREMENTALUPDATE_HH
#include "Update.hh"

class IncrementalUpdate : public Update
{
public:
    IncrementalUpdate();
    ~IncrementalUpdate();
    void pick_dof(int i);
    inline void set_scale(double s) {scl=s;}
    int perform();
    int revert();
private:
    double scl;
    DOF_Info dof;
    int dirn;
};

#endif // INCREMENTALUPDATE_HH
