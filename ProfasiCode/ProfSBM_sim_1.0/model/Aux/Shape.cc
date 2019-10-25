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

#include "Shape.hh"
#include <numeric>

using namespace prf;
using std::vector;

using namespace prf_utils;

void Shape::find_cm()
{
    cm=Vector3();
    meansq=0;
    cm=accumulate(pnts.begin(),pnts.end(),Vector3());
    cm*=(1.0/nref);

    for (int i=0;i<nref;++i) {
        meansq+=(pnts[i]-cm).mag2();
    }

    meansq/=nref;
}

void Shape::find_and_move_to_cm()
{
    meansq=0;
    cm=accumulate(pnts.begin(),pnts.end(),Vector3());
    cm*=(1.0/nref);

    for (int i=0;i<nref;++i) {
        pnts[i]-=cm;
        meansq+=pnts[i].mag2();
    }

    meansq/=nref;
}

Shape::Shape() :nref(0),meansq(0) {}

Shape::~Shape() {}

Shape::Shape(unsigned int i) : nref(i),meansq(0),pnts(i,Vector3()) {}

void Shape::Resize(unsigned int i) {pnts.resize(i,Vector3());nref=i;}

void Shape::Translate(const Vector3 &t)
{
    for (int i=0;i<nref;++i) {
        pnts[i]+=t;
    }
}

void Shape::Rotate(double th, const Vector3 & ax)
{
    for (int i=0;i<nref;++i) {
        Vector3 vv(pnts[i]-cm);
        vv.rotate(th,ax);
        pnts[i]=cm+vv;
    }
}

void Shape::Scale(double xscl,double yscl, double zscl)
{
    for (int i=0;i<nref;++i) {
        pnts[i]=Vector3(xscl*pnts[i].x(),yscl*pnts[i].y(),zscl*pnts[i].z());
    }
}

void Shape::Apply1(Matrix<double> & mtR)
{
    for (int i=0;i<nref;++i) {
        Vector3 vv(pnts[i]);
        pnts[i].x(mtR[0][0]*vv.x()+mtR[0][1]*vv.y()+mtR[0][2]*vv.z());
        pnts[i].y(mtR[1][0]*vv.x()+mtR[1][1]*vv.y()+mtR[1][2]*vv.z());
        pnts[i].z(mtR[2][0]*vv.x()+mtR[2][1]*vv.y()+mtR[2][2]*vv.z());
    }
}

void Shape::Apply2(Matrix<double> & mtR)
{
    for (int i=0;i<nref;++i) {
        Vector3 vv(pnts[i]);
        pnts[i].x(mtR[0][0]*vv.x()+mtR[1][0]*vv.y()+mtR[2][0]*vv.z());
        pnts[i].y(mtR[0][1]*vv.x()+mtR[1][1]*vv.y()+mtR[2][1]*vv.z());
        pnts[i].z(mtR[0][2]*vv.x()+mtR[1][2]*vv.y()+mtR[2][2]*vv.z());
    }
}

void Shape::Save(const char *fln)
{
    Output fout(fln);

    for (int i=0;i<nref;++i) {
        fout<<pnts[i].x()<<"\t"<<pnts[i].y()<<"\t"<<pnts[i].z()<<"\n";
    }

    fout.close();
}

