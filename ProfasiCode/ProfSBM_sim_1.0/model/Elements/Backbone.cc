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

#include "Backbone.hh"
#include "../AA/AminoAcid.hh"

using std::vector;

using namespace UnivConstants;

using namespace prf;

void Backbone::initCommon()
{
    theta[0]=121.7*pi/180;
    theta[1]=111.0*pi/180;
    theta[2]=116.6*pi/180;
    b[0]=1.33;
    b[1]=1.46;
    b[2]=1.52;
    bbf[0].Initialize(b[1],b[2],b[0],pi-theta[1],pi-theta[2]);
    bbf[1].Initialize(b[2],b[0],b[1],pi-theta[2],pi-theta[0]);
    bbf[2].Initialize(b[0],b[1],b[2],pi-theta[0],pi-theta[1]);
    bbb[0].Initialize(b[2],b[1],b[0],pi-theta[1],pi-theta[0]);
    bbb[1].Initialize(b[0],b[2],b[1],pi-theta[2],pi-theta[1]);
    bbb[2].Initialize(b[1],b[0],b[2],pi-theta[0],pi-theta[2]);
}

Backbone::Backbone(int naaa)
{
    initCommon();
    numAminoAcids(naaa);
}

void Backbone::numAminoAcids(int naaa)
{
    naa=naaa;
    nmid=3*naa/2;
    bbatom.resize(3*naa);
    the_bond.resize(3*naa+1,Vector3(0,0,0));
    phi.resize(3*naa,pi);
    bbdof.resize(2*naa,0);
    phitodof.resize(3*naa,-1500);
    frzn.resize(3*naa,false);
    ndof=0;
}

Backbone::Backbone(const Backbone &gbb)
{
    initCommon();
    naa=gbb.naa;ndof=gbb.ndof;nmid=gbb.nmid;
    bbatom=gbb.bbatom;
    the_bond=gbb.the_bond;
    phi=gbb.phi;
    bbdof=gbb.bbdof;
    phitodof=gbb.phitodof;
    frzn=gbb.frzn;
}

void Backbone::registerAA(int i, AminoAcid *aac)
{
    bbatom[3*i]=aac->Nitrogen();
    bbatom[3*i+1]=aac->Calpha();
    bbatom[3*i+2]=aac->Cprime();

    if (aac->OLC()!=P&&aac->OLC()!=DPR) {
        phitodof[3*i]=ndof;
        bbdof[ndof++]=3*i;
        phi[3*i]=-pi;
    } else if (aac->OLC()==P) {
        phi[3*i]=-65.0*pi/180;
        frzn[3*i]=true;
    } else {
        phi[3*i]=+81.0*pi/180;
        frzn[3*i]=true;
    }

    if (!aac->isCis()) phi[3*i+2]=pi; else phi[3*i+2]=0;
    frzn[3*i+2]=true;

    phitodof[3*i+1]=ndof;bbdof[ndof++]=3*i+1;phi[3*i+1]=pi;

    aac->AttachBB(this,i);
}

void Backbone::set_atom_offset(int new_n0_uid)
{
    int old_n0_uid=bbatom[0].UniqueId();
    int offset=new_n0_uid-old_n0_uid;
    if (offset) {
        for (size_t i=0;i<bbatom.size();++i)
            bbatom[i].UniqueId(offset+bbatom[i].UniqueId());
    }
}

Backbone::~Backbone() {}

void Backbone::clear() {
    bbatom.clear();
    the_bond.clear();
    phi.clear();
    bbdof.clear();
    phitodof.clear();
}

int Backbone::AxisAtoms(Atom &a0,Atom &a1,int iloc)
{
    iloc=iloc%ndof;
    a0=bbatom[bbdof[iloc]];a1=bbatom[bbdof[iloc]+1];
    if (bbdof[iloc]<nmid) return 1; else return 0;
}

int Backbone::LocateAA(double xx)
{
    int i=(int)(ndof*xx);
    i=i%ndof;
    i=bbdof[i];
    return i/3;
}

int Backbone::incrDOF(double lx,double xval)
{
    int i=DOFno(lx);
    i=i%ndof;
    i=bbdof[i];
    phi[i]+=xval;

    if (phi[i]<0) phi[i]+=twoPi;

    if (phi[i]>=twoPi) phi[i]-=twoPi;

    return 0;
}

void Backbone::randomize(RandomNumberBase *rangen)
{
    for (int i=0;i<ndof;++i) DOF(i,twoPi*(rangen->shoot()));

    forwardReconstruct(0,3*naa+1);
}

void Backbone::SetToAlphaHelix()
{
    int j=0;

    for (int i=0;i<naa;++i) {
        prf::clog<<"setting amino acid "<<i<<" to helix \n";

        if ((3*i)==bbdof[j]) DOF(j++,(300.0/360.0)*twoPi);

        DOF(j++,(313.0/360.0)*twoPi);

        prf::clog<<"number of dofs fixed = "<<j<<"\n";
    }

    forwardReconstruct(0,3*naa+1);
}

void Backbone::SetToBetaStrand()
{
    int j=0;

    for (int i=0;i<naa;++i) {
        if ((3*i)==bbdof[j]) DOF(j++,(240.0/360.0)*twoPi);

        DOF(j++,(120.0/360.0)*twoPi);
    }

    forwardReconstruct(0,3*naa+1);
}

void Backbone::reconstStart()
{
    Vector3 xcap(1,0,0),zcap(0,0,1);

    if ((bbatom[0].Pos()-bbatom[1].Pos()).mag()>1e-7) {
        xcap=bbatom[1].Pos()-bbatom[0].Pos();
        zcap=bbatom[2].Pos()-bbatom[1].Pos();
        xcap.mag(1);
        zcap-=zcap.dot(xcap)*xcap;
        zcap.mag(1);
    } else bbatom[0].Pos(Vector3(0,0,0));

    the_bond[1]=b[1]*xcap;

    the_bond[2]=b[2]*(-cos(theta[1])*xcap+sin(theta[1])*zcap);

    the_bond[0]=bbb[0].fresh_eval(the_bond[2],the_bond[1],-phi[0]);

    bbatom[1].Pos(bbatom[0].Pos()+the_bond[1]);

    bbatom[2].Pos(bbatom[1].Pos()+the_bond[2]);
}

void Backbone::reconstLast()
{
    int l=3*naa-1;
    Vector3 xcap(1,0,0),zcap(0,0,1);

    if ((bbatom[l].Pos()-bbatom[l-1].Pos()).mag()>1e-7) {
        xcap=bbatom[l].Pos()-bbatom[l-1].Pos();
        zcap=bbatom[l-1].Pos()-bbatom[l-2].Pos();
        xcap.mag(1);
        zcap-=zcap.dot(xcap)*xcap;
        zcap.mag(1);
    } else bbatom[l].Pos(Vector3(0,0,0));

    the_bond[l]=b[2]*xcap;

    the_bond[l-1]=b[1]*(-cos(theta[1])*xcap+sin(theta[1])*zcap);

    the_bond[l+1]=bbf[0].fresh_eval(the_bond[l-1],the_bond[l],phi[l-1]);

    bbatom[l-1].Pos(bbatom[l].Pos()-the_bond[l]);

    bbatom[l-2].Pos(bbatom[l-1].Pos()-the_bond[l-1]);
}

void Backbone::forwardReconstruct(int strt, int iend)
{
    int i=strt;

    if (i<=2) {reconstStart();i=2;}
    else {
        the_bond[i-2]=bbatom[i-2].Pos()-bbatom[i-3].Pos();
        the_bond[i-1]=bbatom[i-1].Pos()-bbatom[i-2].Pos();
    }

    for (;i<iend;++i) {
        the_bond[i]=bbf[i%3].fresh_eval(the_bond[i-2],the_bond[i-1],phi[i-2]);
    }

    for (int j=strt;j<iend;++j) {
        if (j>0 && j<((int)bbatom.size())) {
            bbatom[j].Pos(bbatom[j-1].Pos()+the_bond[j]);
        }
    }
}

void Backbone::reverseReconstruct(int strt, int iend)
{
    int i=strt;

    if (i>=(3*naa-2)) {reconstLast();i=3*naa-2;}
    else {
        the_bond[i+2]=bbatom[i+2].Pos()-bbatom[i+1].Pos();
        the_bond[i+1]=bbatom[i+1].Pos()-bbatom[i].Pos();
    }

    for (;i>=iend;--i) {
        the_bond[i]=bbb[i%3].fresh_eval(the_bond[i+2],the_bond[i+1],-phi[i]);
    }

    for (int j=strt;j>=iend;--j) {
        if (j>0 && j<((int)bbatom.size()))
            bbatom[j-1].Pos(bbatom[j].Pos()-the_bond[j]);
    }
}

void Backbone::FixAngles()
{
    for (int i=0;i<((int)phi.size());++i) {
        while (phi[i]>=pi) phi[i]-=twoPi;

        while (phi[i]<-pi) phi[i]+=twoPi;
    }
}

void Backbone::reconst_bond_vectors()
{
    freconst_bond_vectors(0);
    the_bond[0]=bbb[0].fresh_eval(the_bond[2],the_bond[1],-phi[0]);
}

void Backbone::freconst_bond_vectors(unsigned int i1)
{
    unsigned int i=i1+1;

    for (;i<bbatom.size();++i) {
        the_bond[i]=bbatom[i].Pos()-bbatom[i-1].Pos();
    }

    the_bond[i]=bbf[0].fresh_eval(the_bond[i-2],the_bond[i-1],phi[i-2]);
}

void Backbone::breconst_bond_vectors(unsigned int i1)
{
    unsigned int i=i1;

    for (;i>0;--i) {
        the_bond[i]=bbatom[i].Pos()-bbatom[i-1].Pos();
    }

    the_bond[0]=bbb[0].fresh_eval(the_bond[2],the_bond[1],-phi[0]);
}

double Backbone::find_theta(int i)
{
    if (i==0) {
        prf::clog<<"theta for i=0 is indeterminate. returning 0\n";
        return 0;
    }

    return pi-the_bond[i+1].angle(the_bond[i]);
}

double Backbone::find_phi(int i)
{
    Vector3 v1,v2,v3;
    double cth,sth,phai;
    v1=(the_bond[i].cross(the_bond[i+1]));
    v2=(the_bond[i+1].cross(the_bond[i+2]));
    v1.mag(1);v2.mag(1);
    cth=v1.dot(v2);
    v3=v1.cross(v2);
    sth=v3.mag();

    if (v3.dot(the_bond[i+1])<0) sth=-sth;

    phai=atan2(sth,cth);

    while (phai<0) phai+=pi2;

    while (phai>pi2) phai-=pi2;

    return phai;
}

int Backbone::calc_torsions(std::vector<bool> &specified)
{
    int nfailed=0;
    for (size_t i=0;i<(bbatom.size()-3);++i) {
        if (!frzn[i+1]) {
            if (specified[bbatom[i].UniqueId()] &&
                specified[bbatom[i+1].UniqueId()] &&
                specified[bbatom[i+2].UniqueId()] &&
                specified[bbatom[i+3].UniqueId()]) {
                the_bond[i+1]=bbatom[i+1].Pos()-bbatom[i].Pos();
                the_bond[i+2]=bbatom[i+2].Pos()-bbatom[i+1].Pos();
                the_bond[i+3]=bbatom[i+3].Pos()-bbatom[i+2].Pos();
                phi[i+1]=find_phi(i+1);
            } else ++nfailed;
        }
    }
    return nfailed;
}

void Backbone::CheckProperties()
{
    prf::cout<<"Checking bond properties ...\n";
    reconst_bond_vectors();
    prf::cout<<"bond length 1 ="<<the_bond[1].mag()<<"\texpected value ="
    <<b[1]<<"\n";
    prf::cout<<"----------------------------------------------------------\n";
    prf::cout<<"bond length 2 ="<<the_bond[2].mag()<<"\texpected value ="
    <<b[2]<<"\n";
    prf::cout<<"theta1-2 = "<<find_theta(1)<<"\texpected value ="
    <<theta[1]<<"\n";
    prf::cout<<"----------------------------------------------------------\n";

    for (unsigned int i=3;i<bbatom.size();++i) {
        prf::cout<<"bond length "<<i<<" ="<<the_bond[i].mag()
        <<"\texpected value ="<<b[i%3]<<"\n";
        prf::cout<<"theta "<<i-1<<" = "<<find_theta(i-1)
        <<"\texpected value ="<<theta[(i+2)%3]<<"\n";
        prf::cout<<"phi "<<i-2<<" = "<<find_phi(i-2)
        <<"\texpected value ="<<phi[i-2]<<"\n";
        prf::cout<<"------------------------------------------------------\n";
    }
}


void Backbone::ExportCrds(Shape &sh) const
{
    if (sh.NPoints()!=(int)bbatom.size()) sh.Resize(bbatom.size());

    for (size_t i=0;i<bbatom.size();++i) sh.Point(i,bbatom[i].Position());
}

void Backbone::ExportCrds(Shape &sh,int iaast, int iaand) const
{
    if (sh.NPoints()!=(3*(iaand-iaast))) sh.Resize(3*(iaand-iaast));

    int k=0;

    for (int i=3*iaast;i<(3*iaand);++i,++k) {
        sh.Point(k,bbatom[i].Position());
    }
}

double Backbone::Radius2OfGyration() const
{
    Vector3 cm(0.,0.,0.);

    for (unsigned int i=0;i<bbatom.size();++i) {
        cm+=bbatom[i].Position();
    }

    cm=(1.0/(bbatom.size()))*cm;

    double rdg=0.;

    for (unsigned int i=0;i<bbatom.size();++i) {
        rdg+=(bbatom[i].Position()-cm).mag2();
    }

    return rdg/(bbatom.size());
}

