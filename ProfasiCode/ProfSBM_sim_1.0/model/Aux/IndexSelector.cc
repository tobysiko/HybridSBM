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

#include "IndexSelector.hh"
#include <cmath>
#include "profasi_io.hh"
IndexSelector::IndexSelector() : offset(0) {}

IndexSelector::IndexSelector(const IndexSelector &is) :offset (is.offset)
{
    p=is.p;
    sp=is.sp;
}

IndexSelector &IndexSelector::operator =(const IndexSelector &is)
                                        {
    if (this!=&is) {
        offset=is.offset;
        p=is.p;
        sp=is.sp;
    }
    return *this;
}

IndexSelector::~IndexSelector() {}

void IndexSelector::set_range(size_t i1, size_t i2)
{
    offset=i1;
    p.resize(i2-i1,0);
    sp.resize(p.size()+1,0);
}

size_t IndexSelector::pick(double y)
{
    if (p.size()<2) return offset;
    if (y<0) y=-y;
    if (y>1) y-= (int)y;
    double x1=0,x2=p.size();
    unsigned x=y*x2;
    while (not (sp[x]<=y && sp[x+1]>y)) {
        if (sp[x+1]<=y) x1=x+1;
        if (sp[x]>y) x2=x;
        x=x1+(x2-x1)*(y-sp[x1])/(sp[x2]-sp[x1]);
    }

    return offset+ (size_t) x;
}

void IndexSelector::set_weights(std::vector<double> &wts)
{
    if (wts.size()<p.size()) p.assign(p.size(),0);
    size_t minsize=std::min(wts.size(),p.size());
    for (size_t i=0;i<minsize;++i) p[i]=fabs(wts[i]);
}

void IndexSelector::init()
{
    double wtmx=0;
    for (size_t i=0;i<p.size();++i) {
        if (fabs(p[i])>wtmx) wtmx=p[i];
    }
    if (wtmx==0) {
        prf::cerr<<"IndexSelector> Site weights are all zero. Using flat "
                <<" distribution. \n";
        p.assign(p.size(),1);
    }

    sp[0]=0;
    for (size_t i=0;i<p.size();++i) {
        sp[i+1]=sp[i]+p[i];
    }
    for (size_t i=0;i<sp.size();++i) sp[i]/=sp.back();
}
