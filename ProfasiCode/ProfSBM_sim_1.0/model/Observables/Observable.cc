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

#include "Observable.hh"
#include "../Aux/fileutils.hh"

using namespace prf;

using namespace prf_utils;
using std::string;
using std::vector;

Observable::Observable() : Named(),obsval(0), ntmp(1), myhis(NULL)
{
    requires_population=true;
    fixed_his=false;
    nohis=true;
    gathstat=true;
    hisresume=false;
    userhisrange=false;
    usernbins=false;
    userbinsz=false;
    xbin0=0.01; //if nothing changes it, make histograms of binsize 0.01
    nbins0=100; //default number of bins when bin size is uncertain
    grdtyp=0;
    log_thres=10;
}

Observable::~Observable()
{
    if (myhis) delete(myhis);
}

void Observable::refresh()
{
    obsval=evaluate();
}

double Observable::evaluate()
{
    return obsval;
}

double Observable::delta(Update *u)
{
    double bkp=obsval;
    obsval=evaluate();
    obsdel=obsval-bkp;
    obsval=bkp;
    return obsdel;
}

void Observable::accept(Update *u)
{
    obsval+=obsdel;
    obsdel=0;
}

void Observable::reject(Update *u)
{
    obsdel=0;
}

void Observable::set(string pars)
{
    pars=trim_str(pars);

    if (pars.empty()) return ;

    vector<string> parparts;

    pars=trim_str(pars);

    split_str<vector<string> >(pars,';',parparts);

    for (size_t i=0;i<parparts.size();++i)
        if (!parparts[i].empty()) usrcmd.push_back(parparts[i]);
}

void Observable::set_logger_threshold(int thr)
{
    log_thres=thr;
}

void Observable::rangeEstimate(double &x1,double &x2)
{
    x1=0;x2=1;
    if (!userbinsz) xbin0=0.01;
}

void Observable::write_rtkey(Output &op) {op<<Name()<<"\n";}

void Observable::write_snapshot(Output &op) {op<<"  "<<obsval;}

void Observable::his_range(double xmn, double xmx)
{
    xmin0=xmn;
    xmax0=xmx;
    if (xmin0>xmax0) std::swap(xmin0,xmax0);
    userhisrange=true;
}

void Observable::his_bin_size(double sz)
{
    xbin0=sz;
}

void Observable::his_fill(int i)
{
    if (myhis and gathstat) myhis->put(obsval,i);
}

void Observable::avg_fill(int i)
{
    if (gathstat) avgs[i].put(obsval);
}

void Observable::avg_write(Output &op)
{
    op<<"############## "<< Name() <<" ##############\n";

    for (int j=0;j<ntmp;++j)
        op<<j<<'\t'<<avgs[j].mean() <<'\t'<<avgs[j].stddev() <<"\n";
}

void Observable::his_reset()
{
    if (myhis) myhis->reset();
}

void Observable::avg_reset()
{
    for (int j=0;j<ntmp;++j) avgs[j].reset();
}

void Observable::his_save()
{
    if (myhis and gathstat) {
        if (!fixed_his) myhis->adjust();
        myhis->Export((oprefx+"his_"+Name()).c_str());
    }
}

int Observable::pop_check()
{
    if (p==NULL) {
        prf::cerr<<Name()<<": No population set. Initialization failed. \n";
        return 0;
    } else return 1;
}

/**
* \page opt_Observable Options for all observables
\li \b histogram:  This ensures that the observable keeps a histogram.
Alternative form: make_histogram. Normally this is not needed, unless histograms
have been globally turned off.
\li \b his_range "his_range 25 100" would set the preliminary range of the
histogram for this observable as 25 to 100. This will be changed during the
run if necessary to accommodate the data, unless the "fixed" keyword is also
given at the end : "his_range 25 100 fixed".
\li \b his_bins "his_bins 12" will set the preliminary number of bins in the
histogram. The number of bins is changed when the range is adjusted during a
simulation. If the range is given as "fixed" with the "his_range" option, the
number of bins also remains fixed.
\li \b his_bin_size "his_bin_size 0.5" will set the bin size used for the
histogram to 0.5. Once set, the bin size does not change upon histogram
self adjustments.

\note If \b his_range is present, commands \b his_bins and \b his_bin_size must
come after \b his_range.
\note \b his_bins and \b his_bin_size are mutually exclusive. Only one of them
should be used for one Observable.

\section examples Examples
new_obs %Rg rg2 of_chain 1 ; his_range 4 20 fixed ; his_bins 100 <br>
This sets up a new measurement of the radius of gyration, for the second chain
in the population, with a histogram in the range 4, 20 with 100 bins. The
histogram will maintain this range and number of bins throughout the run,
without trying to adjust it. Out of range points will be lost.<br><br>
new_obs %Rg rg0  of chain 0; his_range 4 20 ; his_bins 100<br>
This sets up a radius of gyration calculation for the first chain with a
histogram that begins collecting data in the range 4 to 20 with 100 bins. But if
the actual data is outside the range, the histogram range and number of bins are
adjusted to cover the data, when that makes sense. <br><br>
new_obs %Rg rg make_histgram <br>
This sets up a radius of gyration of all the atoms present in the population,
and makes a histogram of the data. <br><br>
new_obs %Rg rg his_range 4.3 20 ; his_bin_size 0.25<br>
This makes sure that histograms of the observable named rg are made with a
bin size of 0.25. If necessary, the given range is symmetrically extended so
as to contain an integral number of bins. In this example, the effective
range would be 4.275 to 20.025.
\sa prf_utils::AdaptiveHis, prf::Observable
*/

int Observable::init_obs()
{
    if (requires_population && pop_check()==0) return 0;

    avgs.resize(ntmp);

    vector<string> parts;

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);

        if (parts[0]==string("his_range") && parts.size()>=3) {
            his_range(strtod(parts[1].c_str(),NULL),
                      strtod(parts[2].c_str(),NULL));
            xbin0=(xmax0-xmin0)/nbins0;

            if (parts.size()>3 &&parts[3]=="fixed") fixed_his=true;
        }

        if (parts[0]==string("his_bins") && parts.size()>=2) {
            nbins0=atoi(parts[1].c_str());
            usernbins=true;
            userbinsz=false;
        }

        if (parts[0]==string("his_bin_size") && parts.size()>=2) {
            xbin0=strtod(parts[1].c_str(),NULL);
            usernbins=false;
            userbinsz=true;
        }

        if (parts[0]==string("histogram") || parts[0]==string("make_histogram"))
            nohis=false;

        if (parts[0]==string("hist_resume")) hisresume=true;

        if ((parts[0]==string("gradient_method"))) {
            grdtyp=atoi(parts[1].c_str());
            if (grdtyp!=2) grdtyp=0;
        }

    }

    return 1;
}

void Observable::his_setup()
{
    if (!(myhis or nohis)) myhis=new AdaptiveHis();

    if (myhis) {
        if (!hisresume || myhis->Import((oprefx+"his_"+Name()).c_str())==0) {
            if (!userhisrange) rangeEstimate(xmin0,xmax0);

            myhis->Range(xmin0,xmax0);
            if (usernbins) {
                myhis->Nbins(nbins0);
                Logger(log_thres)<<"Initializing histogram of "<<Name()
                        <<" with range ("<<xmin0<<", "<<xmax0<<"), and "
                        <<nbins0<<" bins.\n";
            } else {
                myhis->set_bin_size(xbin0);
                xmin0=myhis->Xmin();
                xmax0=myhis->Xmax();
                Logger(log_thres)<<"Initializing histogram of "<<Name()
                        <<" with range ("<<xmin0<<", "<<xmax0<<"), and bin size "
                        <<xbin0<<".\n";
            }

            myhis->NBlocks(ntmp);
            myhis->Name(Name());
            myhis->init();
        } else Logger(log_thres)<<"Resumed collecting data for histogram of "
            <<Name()<<"\n";

        if (fixed_his) myhis->disable_adjust();
    }
}

void Observable::set_n_temp(int i)
{
    ntmp=i;
    avgs.resize(ntmp);
}

double Observable::gradientXYZ(std::valarray<double> &ans)
{
    ans=0;
    return 0;
}

double Observable::gradientDOF(std::valarray<double> &ans,
                               std::vector<int> &indxs)
{
    ans=0;
    return 0;
}
