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

#include "OligoOrient.hh"

using namespace prf;

using namespace prf_utils;

using namespace UnivConstants;
using std::string;
using std::vector;

OligoOrient::OligoOrient()
{
    Name("OligoOrient");
    copang=cos(pi/6);
    hbcut=4.6;
    minbt=0.5;
    chbeta.clear();
}

OligoOrient::~OligoOrient()
{
    if (hisb) delete(hisb);
}

/**
* \page opt_OligoOrient OligoOrient: Relative orientation of adjacent chains in Oligomers.
\section options Available options
<ul>
<li><b>hb_cut</b> hb_cut -4.6<br>Minimum interchange hydrogen bond energy for two chains to be considered as parf of one oligomer. </li>
<li><b>beta_cut</b> beta_cut 0.5<br>Minimum beta strand content for a chain to be considered for possible participation in beta sheets. </li>
<li><b>max_angle</b> max_angle 0.5213 <br>Maximum angle between the end to end vectors of two adjacent chains in a beta sheet for them to be stil considered parallel. </li>
<li><b>max_angle_deg</b> max_angle_deg 30 <br> Same as "max_angle", but the values are entered in degrees.  </li>
</ul>

\sa prf::OligoOrient
*/

int OligoOrient::init_obs()
{
    if (Observable::init_obs()==0) return 0;

    if (p->NumberOfChains()<2) {
        prf::cerr<<"The population has only one chain. Oligomer orientation\n"
        <<"observables do not make sense.\n"
        <<Name()<<"> Initialisation failed.\n";
        return 0;
    }

    vector<string> parts;

    //The contact type has to be set explicitly by a user command

    for (size_t i =0;i<usrcmd.size();++i) {
        prf::cout<<"Processed "<<usrcmd[i]<<"\n";
        getchar();
        parts.clear();
        split(usrcmd[i],parts);

        if (parts[0]==string("hb_cut") && parts.size()>=2) {
            minHB(strtod(parts[1].c_str(),NULL));
        }

        if (parts[0]==string("beta_cut") && parts.size()>=2) {
            minBeta(strtod(parts[1].c_str(),NULL));
        }

        if (parts[0]==string("max_angle") && parts.size()>=2) {
            setAngleCut(strtod(parts[1].c_str(),NULL));
        }

        if (parts[0]==string("max_angle_deg") && parts.size()>=2) {
            setAngleCut(strtod(parts[1].c_str(),NULL)*pi/180.0);
        }
    }

    avgb.resize(ntmp);

    chbeta.resize(p->NumberOfChains());

    for (unsigned i=0;i<chbeta.size();++i) {
        char num[100];
        sprintf(num,"%u",i);
        chbeta[i].Name("beta_for_chain_"+string(num));
        chbeta[i].set("default_limits beta");
        chbeta[i].set("nohis");
        sprintf(num,"residues_of_chain %u",i);
        chbeta[i].set(string(num));
        chbeta[i].init_obs();
    }

    cf.set_cutoff(hbcut);

    cf.init(p);
    return 1;
}

void OligoOrient::his_fill(int curT)
{
    if (gathstat) {
        if (myhis) myhis->put(np,curT);
        if (hisb) hisb->put(nm,curT);
    }
}

void OligoOrient::avg_fill(int curT)
{
    if (gathstat) {
        avgs[curT].put(np);
        avgb[curT].put(nm);
    }
}

void OligoOrient::his_reset()
{
    if (myhis) {
        myhis->reset();
        hisb->reset();
    }
}

void OligoOrient::avg_reset()
{
    for (int j=0;j<ntmp;++j) {
        avgs[j].reset();
        avgb[j].reset();
    }
}

void OligoOrient::his_save()
{
    if (gathstat) {
        if (myhis) myhis->Export((oprefx+Name()+string(".PLUS")).c_str());
        if (hisb) hisb->Export((oprefx+Name()+string(".MINUS")).c_str());
    }
}

void OligoOrient::avg_write(Output &op)
{
    op<<"######### "<< Name() <<".PLUS #########\n";

    for (int j=0;j<ntmp;++j)
        op<<j<<'\t'<<avgs[j].mean() <<'\t'<<avgs[j].stddev() <<"\n";

    op<<"######### "<< Name() <<".MINUS ########\n";

    for (int j=0;j<ntmp;++j)
        op<<j<<'\t'<<avgb[j].mean() <<'\t'<<avgb[j].stddev() <<"\n";
}

void OligoOrient::write_rtkey(Output &op)
{
    op<<(Name()+string(".PLUS"))<<"\n";
    op<<(Name()+string(".MINUS"))<<"\n";
}

void OligoOrient::write_snapshot(Output &op)
{
    op<<"  "<<np<<"  "<<nm<<"\n";
}


void OligoOrient::rangeEstimate(double & x1, double & x2)
{
    x1=-0.5;
    x2=p->NumberOfChains()+0.5;
    if (!usernbins) nbins0=p->NumberOfChains()+1;
    if (!userbinsz) xbin0=1;
}

double OligoOrient::evaluate()
{
    np=0,nm=0;
    int NC=p->NumberOfChains();

    for (size_t i=0;i<chbeta.size();++i) chbeta[i].evaluate();

    for (int ich=0;ich<NC;++ich) {
        if (chbeta[ich].Value() <minbt) continue;

        for (int jch=ich+1;jch<NC;++jch) {
            if (chbeta[jch].Value()<minbt) continue;

            if (cf(ich,jch)) {
                Vector3 v1=(p->Chain(ich)->EndToEndVector());
                Vector3 v2=(p->Chain(jch)->EndToEndVector());
                double v1dotv2=(v1.dot(v2))/(v1.mag()*v2.mag());
                if (v1dotv2>copang) ++np; else if (v1dotv2<-copang) ++nm;
            }
        }
    }

    return np;
}

void OligoOrient::his_setup()
{
    if (!userhisrange) rangeEstimate(xmin0,xmax0);

    if (!(myhis or nohis)) myhis=new AdaptiveHis();

    if (myhis) {
        myhis->Range(xmin0,xmax0);
        myhis->Nbins(nbins0);
        myhis->NBlocks(ntmp);
        myhis->Name(Name()+".PLUS");
        myhis->init();
    }

    if (!(hisb or nohis)) hisb=new AdaptiveHis();

    if (hisb) {
        hisb->Range(xmin0,xmax0);
        hisb->Nbins(nbins0);
        hisb->NBlocks(ntmp);
        hisb->Name(Name()+".MINUS");
        hisb->init();
    }
}
