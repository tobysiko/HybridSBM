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

#include "FindCoord.hh"
#include "Constants.hh"

using namespace prf;

using namespace UnivConstants;
FindCoord::FindCoord() {}

FindCoord::FindCoord(double gu,double gv,double gw,double gthuv,double gthvw)
{
    wu=(gw*sin(gthvw))/(gu*sin(gthuv));
    wv1=(gw/gv)*cos(gthvw);
    wv2=-(gw/gv)*sin(gthvw)/tan(gthuv);
    wvxu=gw*sin(gthvw)/(gu*gv*sin(gthuv));
    wsthvw=gw*sin(gthvw);wcthvw=gw*cos(gthvw);
}

FindCoord::FindCoord(double gu,double gv,double gw,double gthuv,
                     double gthvw,double gh)
{
    FindCoord(gu,gv,gw,gthuv,gthvw);
    FixPhi(gh);
}

FindCoord::~FindCoord() {}

void FindCoord::Initialize(double gu,double gv,double gw,double gthuv,double gthvw)
{
    wu=(gw*sin(gthvw))/(gu*sin(gthuv));
    wv1=(gw/gv)*cos(gthvw);
    wv2=-(gw/gv)*sin(gthvw)/tan(gthuv);
    wvxu=gw*sin(gthvw)/(gu*gv*sin(gthuv));
    wsthvw=gw*sin(gthvw);wcthvw=gw*cos(gthvw);
}

void FindCoord::Initialize(double u,double v,double w,double thuv,double thvw,double gh)
{
    Initialize(u,v,w,thuv,thvw);
    FixPhi(gh);
}

Vector3 FindCoord::operator()(const Vector3 &U,const Vector3 &V, double gph)
{
    cph=-cos(gph);sph=-sin(gph);
    return (wu*cph)*U+(wv1+wv2*cph)*V+(wvxu*sph)*(V.cross(U));
}

Vector3 FindCoord::operator()(const Vector3 &U,const Vector3 &V,
                              const Vector3 &VU, double gph)
{
    cph=-cos(gph);sph=-sin(gph);
    return (wu*cph)*U+(wv1+wv2*cph)*V+(wvxu*sph)*VU;
}

void FindCoord::FixPhi(double gph)
{
    sph=sin(gph-pi);cph=cos(gph-pi);
    wu*=cph;wv=wv1+wv2*cph;wvxu*=sph;
}

Vector3 FindCoord::operator()(const Vector3 &U, const Vector3 &V)
{
    return wu*U+wv*V+wvxu*(V.cross(U));
}

Vector3 FindCoord::operator()(const Vector3 &U,const Vector3 &V,
                              const Vector3 &VU)
{
    return wu*U+wv*V+wvxu*VU;
}

Vector3 FindCoord::fresh_eval(Vector3 U,Vector3 V, double gph)
{
    //same calculation as in operator() but without using any stored
    //information about length of u and v and the angle thuv. A lot more
    //forgiving about input vectors U and V. No matter what it gets, so long as
    //they are not colinear, this function will create a vector W with length w
    //(as given during initialization), at an angle thvw to V, and with a torsion
    //gph around V relative to the direction of U.
    cph=-cos(gph);sph=-sin(gph);
    U.mag(1);V.mag(1);
    Vector3 UcrossV=U.cross(V);
    double udotv=fabs(U.dot(V)),ucrossv=UcrossV.mag();
    return (wsthvw*cph/ucrossv)*U
           +(wcthvw-wsthvw*udotv*cph/ucrossv)*V
           -(wsthvw*sph/ucrossv)*(UcrossV);
}
