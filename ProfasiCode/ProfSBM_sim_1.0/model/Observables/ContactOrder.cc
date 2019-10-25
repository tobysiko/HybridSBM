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

#include "ContactOrder.hh"

using namespace prf_utils;
using std::string;
using std::vector;
ContactOrder::ContactOrder()
{
    Name("ContactOrder");
    ich=0;
    dcut=6.0;
    cf.min_seq_sep(1);
    cf.set_cutoff(dcut);
//    cf.minimum_links(1);
}

ContactOrder::~ContactOrder() {}

/**
\page opt_ContactOrder Relative Contact Order
\section options Available options
<ul>
<li><b>of_chain</b> of_chain 3<br> Associate with one particular chain in a multi-chain simultion.</li>
<li><b>distance_cut</b> distance_cut 6.5<br> Distance cut off used to decide whether two residues are in contact in a contact order calculation. Default is 6.0 Angstroems.</li>
<li><b>min_seq_sep</b> min_seq_sep 2<br> Minimum sequence separation between two residues for them to be considered for contact order calculations. The default is 1, meaning even adjacent residues will contribute. </li>
</ul>
\section examples Example
new_obs %ContactOrder CO4 of_chain 3 ; distance_cut 6.5 ;
\sa ContactOrder
*/
int ContactOrder::init_obs()
{
    Logger blog(log_thres);

    if (Observable::init_obs()==0) return 0;

    vector<string> parts;

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);

        if (parts[0]==string("of_chain") && parts.size()>=2) {
            int tmpint=atoi(parts[1].c_str());

            if ((tmpint<0) || (tmpint > p->NumberOfChains())) tmpint=0;

            of_chain(tmpint);

            blog<<"ContactOrder("<<Name()<<")> associating with chain "
            <<tmpint<<"\n";
        }

        if (parts[0]==string("distance_cut") && parts.size()>=2) {
            set_cutoff(strtod(parts[1].c_str(),NULL));
            blog<<"ContactOrder("<<Name()<<")> Set distance cut to "
            <<dcut<<"\n";
            cf.set_cutoff(dcut);
        }

        if (parts[0]==string("min_seq_sep") && parts.size()>=2) {
            cf.min_seq_sep(atoi(parts[1].c_str()));
            blog<<"ContactOrder("<<Name()<<")> Set minimum sequence "
            <<"separation to "<<atoi(parts[1].c_str())<<"\n";
        }
    }

    cf.init(p);

    return 1;
}

double ContactOrder::evaluate()
{
    int nr1=0,nr2=p->Chain(ich)->numLigands();

    int nterms=0;
    double dsavg=0;

    for (int i=nr1;i<nr2;++i) {
        for (int j=i+1;j<nr2;++j) {
            if (cf(i,j)) {
//                prf::cout<<"Contributing difference "<<abs(i-j)<<"\n";
                dsavg+=abs(i-j);
                ++nterms;
            }
        }
    }

    if (nterms) dsavg/=nterms;

//    prf::cout<<"Number of contributing terms = "<<nterms<<"\n";
    return 100.0*dsavg/(nr2-nr1);
}

void ContactOrder::rangeEstimate(double &x1,double &x2)
{
    x1=0;
    x2=30.0;
    if (!userbinsz) xbin0=0.5;
}
