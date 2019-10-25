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

#include "rmsd.hh"
#include "profasi_io.hh"
#include "Constants.hh"

using namespace prf;

using namespace prf_utils;

using namespace UnivConstants;
using std::max;
using std::min;

RMSD::RMSD() {}

RMSD::~RMSD() {}

double RMSD::operator()(Shape &sh1, Shape &sh2, bool sh1incms)
{
    Vector3 v1,v2;
    int npt=min(sh1.NPoints(),sh2.NPoints());
    npt=max(1,npt);
    if (!sh1incms) sh1.find_and_move_to_cm();
    sh2.find_cm();
    double gyr=(sh1.gyr_r2()+sh2.gyr_r2());

    std::valarray<double> AA(9);

    for (int k=0;k<npt;++k) {
        v1=sh1.Point(k);v2=sh2.Point(k);
        AA[0]+=v1.x()*v2.x();
        AA[1]+=v1.x()*v2.y();
        AA[2]+=v1.x()*v2.z();
        AA[3]+=v1.y()*v2.x();
        AA[4]+=v1.y()*v2.y();
        AA[5]+=v1.y()*v2.z();
        AA[6]+=v1.z()*v2.x();
        AA[7]+=v1.z()*v2.y();
        AA[8]+=v1.z()*v2.z();
    }

    AA/=npt;

    f.set_position(AA);
    return sqrt(fabs(gyr-2*f.value()));
}

double RMSD::eval_w_grad(Shape &sh1, Shape &sh2, std::valarray<double> &dr)
{
    Vector3 v1,v2;
    int npt=min(sh1.NPoints(),sh2.NPoints());
    npt=max(1,npt);
    dr=0;
    sh1.find_and_move_to_cm();
    sh2.find_and_move_to_cm();
    double gyr=(sh1.gyr_r2()+sh2.gyr_r2());

    Matrix<double> dA(3*npt,9);
    std::valarray<double> AA(9);

    for (int k=0;k<npt;++k) {
        v1=sh1.Point(k);v2=sh2.Point(k);
        AA[0]+=v1.x()*v2.x();
        AA[1]+=v1.x()*v2.y();
        AA[2]+=v1.x()*v2.z();
        AA[3]+=v1.y()*v2.x();
        AA[4]+=v1.y()*v2.y();
        AA[5]+=v1.y()*v2.z();
        AA[6]+=v1.z()*v2.x();
        AA[7]+=v1.z()*v2.y();
        AA[8]+=v1.z()*v2.z();
        for (int j=0;j<3;++j) {
            for (int i=0;i<3;++i) {
                dA.set(3*k+i,3*j+i,v1[j]/npt);
            }
            dr[3*k+j]=2*v2[j]/npt;
        }
    }

    AA/=npt;

    f.set_position(AA);

    double ans=sqrt(fabs(gyr-2*f.value()));

    std::valarray<double> dL(0.0,9);
    f.gradient(dL);

    for (int i=0;i<3*npt;++i) {
        for (int j=0;j<9;++j) {
            dr[i]-=2*dA.get(i,j)*dL[j];
        }
        dr[i]/=(2*ans);
    }
    return ans;
}

rmsd_help_fn::rmsd_help_fn() : NDFunctional()
{
    nd=9;
    cur.resize(9,0);
}

rmsd_help_fn::~rmsd_help_fn() {}

double rmsd_help_fn::value()
{
    double a,b,c,q,r,theta,w[3];

    double detA=cur[0]*(cur[4]*cur[8]-cur[7]*cur[5])
                -cur[1]*(cur[3]*cur[8]-cur[5]*cur[6])
                +cur[2]*(cur[3]*cur[7]-cur[6]*cur[4]);
    double AtA[9]={0,0,0,0,0,0,0,0,0};

    for (int i=0;i<3;++i) {
        for (int k=0;k<3;++k) {
            for (int j=0;j<3;++j) {
                AtA[3*i+j]+= (cur[3*k+i]*cur[3*k+j]);
            }
        }
    }

    a=-(AtA[0]+AtA[4]+AtA[8]);

    b=AtA[0]*AtA[4]+AtA[0]*AtA[8]+AtA[4]*AtA[8]-
      AtA[1]*AtA[1]-AtA[2]*AtA[2]-AtA[5]*AtA[5];
    c=-detA*detA;
    q=(a*a-3*b)/9;
    r=(2*a*a*a-9*a*b+27*c)/54;

    if ((r*r/(q*q*q))>=1) {
        prf::cerr<<"rmsd: error r^2="<<r*r<<" q^3="
        <<q*q*q<<"ratio="<<r*r/(q*q*q)<<"\n";
        return 0;
    }

    q=sqrt(q);

    theta=acos(r/(q*q*q));

    w[0]=sqrt(fabs(-2*q*cos(theta/3)-a/3));
    w[1]=sqrt(fabs(-2*q*cos((theta+pi2)/3)-a/3));
    w[2]=sqrt(fabs(-2*q*cos((theta-pi2)/3)-a/3));

    if (detA<0) {
        int tmp = w[1] > w[0] ? 1 : 0;
        return (w[tmp]+fabs(w[(tmp+1)%3]-w[(tmp+2)%3]));
    } else return (w[0]+w[1]+w[2]);
}
