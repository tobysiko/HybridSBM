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

#include "NDFunctional.hh"
#include "LineProjection.hh"
#include "profasi_io.hh"
#include <cmath>


NDFunctional::NDFunctional() : nd(0), h_small(0.00001) {}
NDFunctional::~NDFunctional() {}

double NDFunctional::operator ()(std::valarray<double> &x)
{
    set_position(x);
    return value();
}

double NDFunctional::value()
{
    return 0;
}

void NDFunctional::allocate_coords(int ndim)
{
    if (ndim>0) {
        cur.resize(ndim,0);
        nd=ndim;
    }
}

void NDFunctional::set_position(std::valarray<double> &x)
{
    cur=x;
}

void NDFunctional::get_position(std::valarray<double> &x)
{
    x=cur;
}

double NDFunctional::get_coordinate(int i)
{
    return cur[i];
}

void NDFunctional::set_coordinate(int i, double x)
{
    cur[i]=x;
}

double NDFunctional::delta(int idir, double h)
{
    double x0=get_coordinate(idir);
    set_coordinate(idir,x0+h);
    double vp,vm;
    vp=value();
    set_coordinate(idir,x0-h);
    vm=value();
    set_coordinate(idir,x0);
    return vp-vm;
}

void NDFunctional::gradient(std::valarray<double> &g)
{
    for (int i=0;i<nd;++i) {
        g[i]=partial_derivative(i);
    }
}

double NDFunctional::partial_derivative(int idir)
{
    double con=1.4, big=1.0e30, safe=2.0;
    double con2=con*con;
    const int ntab=10;
    double err,errt,e1,e2,fac,hh,ans=0;
    double a[ntab][ntab];
    hh=h_small;
    a[0][0]=delta(idir,hh)/(2*hh);
    err=big;
    for (int i=1;i<ntab;++i) {
        hh/=con;
        a[i][0]=delta(idir,hh)/(2*hh);
        fac=con2;
        for (int j=1;j<=i;++j) {
            a[i][j]=(a[i][j-1]*fac-a[i-1][j-1])/(fac-1);
            fac=con2*fac;
            e1=fabs(a[i][j]-a[i][j-1]);
            e2=fabs(a[i][j]-a[i-1][j-1]);
            errt=(e1>e2)?e1:e2;
            if (errt<=err) {
                err=errt;
                ans=a[i][j];
            }
        }
        if (fabs(a[i][i]-a[i-1][i-1])>=safe*err) break;
    }
    return ans;
}
