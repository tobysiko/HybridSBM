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

#include "LineProjection.hh"

LineProjection::LineProjection() {}

LineProjection::~LineProjection() {}

double LineProjection::operator ()(double t)
{
    fpos=orig+t*dirc;
    return (*ndf)(fpos);
}

void LineProjection::set_base_function(NDFunctional *gf)
{
    if (gf!=NULL) {
        ndf=gf;
        orig.resize(gf->n_dim(),0);
        fpos.resize(gf->n_dim(),0);
        dirc.resize(gf->n_dim(),0);
        dirc[0]=1;
    }
}

void LineProjection::set_origin(std::valarray<double> &org)
{
    size_t mindm=org.size();
    if (orig.size()<mindm) mindm=orig.size();
    for (size_t i=0;i<mindm;++i) orig[i]=org[i];
}

void LineProjection::set_direction(std::valarray<double> &drc)
{
    size_t mindm=drc.size();
    if (dirc.size()<mindm) mindm=dirc.size();
    double vmag=0;
    for (size_t i=0;i<dirc.size();++i) {
        if (i<mindm) dirc[i]=drc[i];
        else dirc[i]=0;
        vmag+=dirc[i]*dirc[i];
    }
    if (vmag!=0) dirc/=sqrt(vmag);
}
