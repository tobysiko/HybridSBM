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

#include "AtomCoordinates.hh"
#include "../Aux/profasi_io.hh"

using std::valarray;

using namespace prf;

using namespace prf_utils;

double AtomCoordinates::BOXL=100;
double AtomCoordinates::HALFBX;
std::valarray<double> AtomCoordinates::atc;
std::valarray<double> AtomCoordinates::btc;
std::valarray<double> * AtomCoordinates::uc=&atc;
int AtomCoordinates::ntotatoms=0;

AtomCoordinates::AtomCoordinates(int i) {Id(i);}

AtomCoordinates::~AtomCoordinates() {}

AtomCoordinates::AtomCoordinates(const AtomCoordinates &cr) {locn=cr.locn;}

void AtomCoordinates::Initialize(int tot_atoms)
{
    HALFBX=0.5*BOXL;
    ntotatoms=tot_atoms;
    atc.resize(tot_atoms*3);
    btc.resize(tot_atoms*3);
    switch_to_regular();
}

void AtomCoordinates::update(int strt, int iend)
{
    strt*=3;
    iend*=3;

    for (int i=strt;i<iend;++i) {
        btc[i]=atc[i];
    }
}

void AtomCoordinates::revert(int strt, int iend)
{
    strt*=3;
    iend*=3;

    for (int i=strt;i<iend;++i) {
        atc[i]=btc[i];
    }
}

Vector3 AtomCoordinates::CenterOfMass(int strt, int iend)
{
    Vector3 ans(0,0,0);

    if (strt>=iend) {
        prf::cerr<<"invalid range of atoms to find center of mass: "
        <<strt<<" to "<<iend<<"\n";
        return ans;
    }

    for (int i=strt;i<iend;++i) ans+=vec(i);

    return (1.0/ (iend-strt)) *ans;
}

Vector3 AtomCoordinates::CenterOfMass3(int i1, int i2, int i3)
{
    Vector3 ans(0,0,0);
    ans = vec(i1) + vec(i2) + vec(i3);

    return 0.25 *ans;
}

void AtomCoordinates::center()
{
    Vector3 cms=CenterOfMass(0,numberOfAtoms());
    BlockTranslate(-cms,0,numberOfAtoms());
}


double AtomCoordinates::rgyr(int istrt, int iend)
{
    if (istrt>=iend) {
        prf::cerr<<"invalid range of atoms to find center of mass: "
        <<istrt<<" to "<<iend<<"\n";
        return 0;
    }

    double vsq=0;

    Vector3 vc;

    for (int i=istrt;i<iend;++i) {
        vc+=vec(i);
        vsq+=vec(i).mag2();
    }

    return (vsq / (iend-istrt)) - ((1.0/(iend-istrt))*vc).mag2();
}

void AtomCoordinates::BlockRotate(double tht, int strt, int iend,
                                  AtomCoordinates ax1, AtomCoordinates ax2)
{
    BlockRotate(tht,strt,iend,ax1.value(),ax2-ax1);
}

void AtomCoordinates::BlockRotate(double tht, int strt, int iend,
                                  Vector3 org, Vector3 axs)
{
    double cdph(cos(tht)),sdph(sin(tht));
    Vector3 un_ax= (1/axs.mag()) *axs;
    Vector3 Ax(un_ax.x() *un_ax);
    Vector3 Ay(un_ax.y() *un_ax);
    Vector3 Az(un_ax.z() *un_ax);
    Vector3 Bx(cdph,-un_ax.z() *sdph,un_ax.y() *sdph);
    Vector3 By(-Bx.y(),cdph,-un_ax.x() *sdph);
    Vector3 Bz(-Bx.z(),-By.z(),cdph);
    Vector3 Cx=Bx+ (1-Bx.x()) *Ax-Bx.y() *Ay-Bx.z() *Az;
    Vector3 Cy=By-By.x() *Ax+ (1-By.y()) *Ay-By.z() *Az;
    Vector3 Cz=Bz-Bz.x() *Ax-Bz.y() *Ay+ (1-Bz.z()) *Az;

    for (int i=strt;i<iend;++i) {
        Vector3 dv=vec(i)-org;
        vec(i,org.x()+Cx.dot(dv),org.y()+Cy.dot(dv),org.z()+Cz.dot(dv));
    }
}

void AtomCoordinates::BlockTranslate(const Vector3 & amnt, int strt, int iend)
{
    strt*=3;iend*=3;

    for (int i=strt;i<iend;) {
        (*uc)[i++]+=amnt.x();
        (*uc)[i++]+=amnt.y();
        (*uc)[i++]+=amnt.z();
    }
}

void AtomCoordinates::EnforceBC(int i1, int i2)
{
    int j;
    double cxmin,cxmax,cymin,cymax,czmin,czmax,xtr=0,ytr=0,ztr=0;
    int minfloor,maxceil;
    BoundingBox(i1,i2,cxmin,cxmax,cymin,cymax,czmin,czmax);

    if (fabs(cxmax-cxmin) >BOXL||
        fabs(cymax-cymin) >BOXL||
        fabs(czmax-czmin) >BOXL) {
        prf::cerr<<"A chain has spread out over a distance larger than a box\n"
        <<"length. Updates on this chain are no longer reliable.\n"
        <<"The relevant chain consists of atoms between unique ids "<<i1<<" and "<<i2<<"\n"
        <<"Perhaps a larger box length should be used. \n"
        <<"current box size = "<<BOXL<<" whereas, chain span is ("
        << (cxmax-cxmin) <<", "<< (cymax-cymin) <<", "<< (czmax-czmin)
        <<"\n";
    }

    minfloor= (int)(floor((cxmin+HALFBX) /BOXL));

    maxceil= (int)(ceil((cxmax-HALFBX) /BOXL));

    if (cxmax>HALFBX) xtr=minfloor*BOXL;

    if (cxmin<-HALFBX) xtr=maxceil*BOXL;

    minfloor= (int)(floor((cymin+HALFBX) /BOXL));

    maxceil= (int)(ceil((cymax-HALFBX) /BOXL));

    if (cymax>HALFBX) ytr=minfloor*BOXL;

    if (cymin<-HALFBX) ytr=maxceil*BOXL;

    minfloor= (int)(floor((czmin+HALFBX) /BOXL));

    maxceil= (int)(ceil((czmax-HALFBX) /BOXL));

    if (czmax>HALFBX) ztr=minfloor*BOXL;

    if (czmin<-HALFBX) ztr=maxceil*BOXL;

    if (fabs(xtr) >0.01 || fabs(ytr)>0.01 || fabs(ztr)>0.01) {      //one if
        i1*=3;
        i2*=3;

        for (j=i1;j<i2;j+=3) {
            atc[j]-=xtr;
            atc[j+1]-=ytr;
            atc[j+2]-=ztr;
        }
    }

}

double AtomCoordinates::total_error()
{
    double disdiff=0;

    for (int i=0;i<numberOfAtoms();++i) {
        disdiff+= (diff_from_backup(i).mag());
    }

    return disdiff;
}

AtomCoordinates & AtomCoordinates::operator= (const AtomCoordinates &g)
{
    if (this!=&g) {
        value(g.x(),g.y(),g.z());
    }

    return *this;
}

Vector3 operator+ (const Vector3 &v3, const AtomCoordinates &g)
{
    return (v3+g.value());
}

Vector3 operator- (const Vector3 &v3, const AtomCoordinates &g)
{
    return (v3-g.value());
}

double AtomCoordinates::per_sep(const double x1, const double x2)
{
    double xdf=x1-x2;

    while (xdf<=-HALFBX) xdf+=BOXL;

    while (xdf>HALFBX) xdf-=BOXL;

    return xdf;

    /*   double ans=xdf-(floor(xdf/BOXLX))*BOXLX; */
    /*   return (ans>HALFBX)?(ans-BOXLX):ans; */
}

double AtomCoordinates::pos_per_sep(const double x1, const double x2)
{
    double xdf=x1-x2;

    while (xdf<0) xdf+=BOXL;

    while (xdf>=BOXL) xdf-=BOXL;

    return xdf;

    /*   double ans=xdf-(floor(xdf/BOXLX))*BOXLX; */
    /*   return (ans>HALFBX)?(ans-BOXLX):ans; */
}

void AtomCoordinates::CheckDiff()
{
    double disdiffavg=0,disdiffmax=0,disdiffcur=0;
    disdiffavg=disdiffmax=0;

    for (int i=0;i<numberOfAtoms();++i) {
        disdiffavg+= (disdiffcur=diff_from_backup(i).mag());

        if (disdiffcur>disdiffmax) disdiffmax=disdiffcur;
    }

    disdiffavg/=numberOfAtoms();

    if (disdiffavg>1e-5 || disdiffmax>1e-6) {
        prf::cerr <<"After reconstruction avg error = "<<disdiffavg
        <<" max error = "<<disdiffmax<<"\n";
    }
}

void AtomCoordinates::write_pdb_box_info(FILE *fp)
{
    double bxlx,bxly,bxlz,scx,scy,scz,thxy=90,thyz=90,thzx=90;
    bxlx=bxly=bxlz=BOXL;
    scx=1/bxlx;
    scy=1/bxly;
    scz=1/bxlz;
    fprintf(fp,"CRYST1  %7.3f  %7.3f  %7.3f %6.2f %6.2f %6.2f P 1    \n"
            ,bxlx,bxly,bxlz,thxy,thyz,thzx);
    fprintf(fp,"ORIGX1      1.000000  0.000000  0.000000        0.00000  \n");
    fprintf(fp,"ORIGX2      0.000000  1.000000  0.000000        0.00000  \n");
    fprintf(fp,"ORIGX3      0.000000  0.000000  1.000000        0.00000  \n");
    fprintf(fp,"SCALE1      %8.6f  0.000000  0.000000        0.50000\n",scx);
    fprintf(fp,"SCALE2      0.000000  %8.6f  0.000000        0.50000\n",scy);
    fprintf(fp,"SCALE3      0.000000  0.000000  %8.6f        0.50000\n",scz);
}

void AtomCoordinates::BoundingBox(int i1, int i2,
                                  double &xmn,
                                  double &xmx,
                                  double &ymn,
                                  double &ymx,
                                  double &zmn,
                                  double &zmx)
{
    xmn=ymn=zmn=1e70;
    xmx=ymx=zmx=-1e70;
    i1*=3;
    i2*=3;

    for (int k=i1;k<i2;++k) {
        if (atc[k]<xmn) xmn=atc[k];

        if (atc[k]>xmx) xmx=atc[k];

        if (atc[++k]<ymn) ymn=atc[k];

        if (atc[k]>ymx) ymx=atc[k];

        if (atc[++k]<zmn) zmn=atc[k];

        if (atc[k]>zmx) zmx=atc[k];
    }
}



