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

#include "UnaryFunction.hh"
#include <cmath>

UnaryFunction::UnaryFunction() : h_small(0.0001) {}

UnaryFunction::~UnaryFunction() {}

double UnaryFunction::operator()(double x)
{
    return x*x;
}

void UnaryFunction::bracket_minimum(double xstart, double step,
                                    double &x0, double &x1, double &x2)
{
    x0=xstart-step;
    x1=xstart;
    double v0,v1,v2;
    v0=(*this)(x0);
    v1=(*this)(x1);
    if (v1>v0) {
        v2=v0;v0=v1;v1=v2;
        x2=x0;x0=x1;x1=x2;
        step=-step;
    }

    do {
        step*=2;
        x2=x1+step;
        v2=(*this)(x2);
        if (v2<=v1) {
            x0=x1;
            x1=x2;
            v0=v1;
            v1=v2;
        } else break;
    } while (1);
}

double UnaryFunction::nearest_minimum(double &x0, double stp)
{
    double x1=0,x2=0;
    bracket_minimum(x0,stp,x0,x1,x2);
    const int itmax=100;
    const double zeps=1.0e-10,cgold=0.3819660;
    const double tol=1.0e-7;

    double a, b,d=0.5,e,etemp, fu, fv, fw, fx, p,q,r,tol1, tol2, u, v,w,x,xm;
    a=(x0<=x2?x0:x2);
    b=(x0>x2?x0:x2);
    w=x=v=x1;
    e=0;
    fv=fw=fx=(*this)(x);

    for (int iter=0; iter<itmax;++iter) {
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+zeps;
        tol2=2.0*tol1;
        if (fabs(x-xm)<=(tol2-0.5*(b-a))) break;
        if (fabs(e)>tol1) {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2*(q-r);
            if (q>0) p=-p;
            q=fabs(q);
            etemp=e;
            e=d;
            if ((fabs(p)>=fabs(0.5*q*etemp) or p<q*(a-x) or p>=q*(b-x))) {
                if (x>=xm) e=a-x; else e=b-x;
                d=cgold*e;
            } else {
                d=p/q;
                u=x+d;
                if (u-a<tol2 or b-u <tol2) d=(xm-x>=0?fabs(tol1):-fabs(tol1));
            }
        } else {
            if (x>=xm) e=a-x; else e=b-x;
            d=cgold*e;
        }

        if (fabs(d)>=tol1) u=x+d; else u=x+(d>=0?fabs(tol1):-fabs(tol1));
        fu=(*this)(u);
        if (fu<=fx) {
            if (u>=x) a=x; else b=x;
            v=w; fv=fw;
            w=x; fw=fx;
            x=u; fx=fu;
        } else {
            if (u<x) a=u; else b=u;
            if (fu<=fw or w==x) {
                v=w; fv=fw;
                w=u; fw=fu;
            } else if (fu<=fv or v==x or v==w) {
                v=u;
                fv=fu;
            }
        }
    }
    x0=x;
    return fx;
}

double UnaryFunction::delta(double x, double h)
{
    return (f(x+h)-f(x-h));
}

double UnaryFunction::derivative(double x)
{
    return derivative(x,h_small);
}

double UnaryFunction::derivative(double x, double h)
{
    double err=0;
    return derivative(x,h,err);
}

double UnaryFunction::derivative(double x, double h, double &err)
{
    double con=1.4, big=1.0e30, safe=2.0;
    double con2=con*con;
    const int ntab=10;
    double errt,e1,e2,fac,hh,ans=0;
    double a[ntab][ntab];
    hh=h;
    a[0][0]=delta(x,hh)/(2*hh);
    err=big;
    for (int i=1;i<ntab;++i) {
        hh/=con;
        a[i][0]=delta(x,hh)/(2*hh);
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

