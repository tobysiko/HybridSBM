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

#include "WangLandau.hh"
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace prf;

using namespace UnivConstants;

WangLandau::WangLandau() : GMC()
{
    converr=0.2;
    gfact0=2.718281828459045;
    init_g_factor();
    nbins=50;
    bt=1.0;
    gvalsread=usr_energy_range=usr_n_bins=using_htmc=false;
    output_prefix="tmp_";
    par.option("energy_range","er",2,"(range of energies for initial histogram)");
    par.option("number_of_bins","nb",1,"(number of bins in initial histogram)");
    par.new_switch("use_htmc","htmc",false,
                   "(Use high temperature MC when energy exceeds upper limit)");
    par.option("flatness_cutoff","fc",1,
               "(noise level below which a histogram is considered flat)");
    par.option("initial_wl_factor","gf0",1,
               "(multiplicative factor for first round of iterations)");
    par.option("read_g_from","g",1,"(read in g values, from saved file)");
    par.disable("max_temperature");
    par.disable("min_temperature");
    par.disable("num_temperatures");
    par.disable("temperature_file");
    descr="Wang-Landau iterations";
}
/**
\page wlopts Extra options for Wang Landau iterations

\li \b --n_stages or \b -ns : How often the Wang Landau multiplicative factor
is updated. Example: --n_stages 8
\li \b --energy_range or \b -er : Initial range for energy histogram in ProFASi
energy units. Example: -er 20 80
\li \b --number_of_bins or \b -nb : Number of bins in the energy histogram at
initialisation. Example: -nb 100
\li \b --flatness_cutoff or \b -fc : The fraction of the population of the most
occupied bin, that the least occupied bin can trail by, for the energy histogram
to be considered flat. Example: -fc 0.2
\li \b --initial_wl_factor or \b -gf0 : The initial value for the multiplicative
factor used to update the density of states in one bin when the simulation
enters it. The default value is e (2.718281828459045) as in the original
Wang Landau article. But it can be set to something else using this option.
\li \b --use_htmc or \b -htmc : Whether or not to use the high temperature
Monte Carlo to decorrlate when the simulation goes beyond the upper bound of
energy. It is not set by default. So, if the energy exceeds the upper bound
the update is rejected, and the system stays where it was, unless you use
this switch. Example : --use_htmc
*/
WangLandau::~WangLandau() {}

int WangLandau::bin_of(double x)
{
    int ibin=floor((x-emin)/ebin);

    if (ibin<0) ibin=0;

    if (ibin>=nbins) ibin=nbins-1;

    return ibin;
}

void WangLandau::init_g()
{
    gvals.clear();
    gvals.resize(nbins,0); // We are actually storing ln(g) values
}

void WangLandau::init_his()
{
    his.clear();
    his.resize(nbins,0);
    ebin=(emax-emin)/nbins;
}

void WangLandau::reset_his()
{
    double mingval=1e100;
    for (int i=0;i<nbins;++i) {
        mingval=std::min(gvals[i],mingval);
        his[i]=0;
    }
    for (int i=0;i<nbins;++i) gvals[i]-=mingval;
}

int WangLandau::extend_till(double lowe)
{
    int nextrabins=0;
    if (lowe<emin) {
        nextrabins=ceil((emin-lowe)/ebin);
        emin-=nextrabins*ebin;
        std::vector<double> hisb,gvb;
        hisb=his;
        gvb=gvals;
        prf::cout<<"While extending to energy "<<lowe<<" adding "<<nextrabins
                <<" bins.\n";
        nbins+=nextrabins;
        init_his();
        init_g();
        for (int i=nextrabins;i<nbins;++i) {
            his[i]=hisb[i-nextrabins];
            gvals[i]=gvb[i-nextrabins];
        }
    }
    return nextrabins;
}

bool WangLandau::converged()
{
    double levellow=1e50,levelhigh=0;
    int npopbins=0;
    for (int i=0;i<nbins;++i) {
        if (his[i]==0) continue; else ++npopbins;
        levellow=std::min(his[i],levellow);
        levelhigh=std::max(his[i],levelhigh);
    }
//    prf::cout<<"Convergence check: range = "<<levellow<<" to "<<levelhigh<<" with "<<npopbins<<" populated bins\n";
    return (npopbins>(0.90*nbins))&&(levellow>=((1-converr)*levelhigh));
}

void WangLandau::init()
{
    if (!gvalsread) {
        if (not usr_energy_range) {
            Etot etot;
            etot.connect(ffh->interaction_potential());
            etot.rangeEstimate(emin,emax);
            prf::cout<<"E range "<<emin<<", "<<emax<<"\n";
            emax+=(2*(emax-emin));
        }
        if (not usr_n_bins) n_bins(40);
        init_g();
        init_his();
    }
    prevbins=nbins;
}

int WangLandau::bring_to_range()
{
    int niter=0;
    double etot=ffh->interaction_potential()->value();
//    prf::cout<<"Bring to range "<<etot<<"\n";
    while (etot>=emax) {
        MC::Step();
        if (niter++%cyclgt==0) ffh->interaction_potential()->reset_total();
        etot=ffh->interaction_potential()->value();
    }

    if (etot<emin) extend_till(etot);
//    prf::cout<<"Bring to range: took "<<niter<<" iterations. Final e = "<<etot<<"\n";
    return niter;
}

int WangLandau::Step()
{
    double enew=0, de=0;
    int ibin=0,jbin=0;

    Update *up=perform_update();

    de=ffh->interaction_potential()->deltaE(up);
    if (de<-100) {
        prf::cerr<<"Large negative energy change for update "<<up->Name()<<" de = "<<de<<"\n";
        for (size_t j=0;j<ffh->interaction_potential()->n_terms();++j) {
            prf::cout<<ffh->interaction_potential()->term(j)->Name()<<": "
                    <<ffh->interaction_potential()->term(j)->deltaE()<<"\n";
        }
    }
    double etot=ffh->interaction_potential()->value();
    enew=etot+de;
//    prf::cout<<"etot = "<<etot<<", de = "<<de<<"\n";
    ibin=bin_of(etot);
    double wt=up->intrinsic_weight();
    if (enew<emax) {
        if (enew<emin) extend_till(enew);
        jbin=bin_of(enew);
        double dg=gvals[ibin]-gvals[jbin];
        double r=ran->shoot();
        if (r<(wt*exp(dg))) {
          //prf::cout<<"Accepting update, method 1\n";
            etot=enew;
            ++his[jbin];
            gvals[jbin]+=lngf;
            ffh->interaction_potential()->accept(up);
            up->accept();
//            getchar();
            return 1;
        }
    } else if (using_htmc && (ran->shoot()<(wt*exp(-bt*(de))))) {
//             prf::cout<<"Accepting update, method 2\n";
        etot=enew;
        ffh->interaction_potential()->accept(up);
        up->accept();
        bring_to_range();
        ibin=bin_of(etot);
        ++his[ibin];
        gvals[ibin]+=lngf;
//        getchar();
        return 1;
    }
    ++his[ibin];
    gvals[ibin]+=lngf;
//     prf::cout<<"Rejecting update \n";
    ffh->interaction_potential()->reject(up);
    up->revert();
//    getchar();
    return 0;
}

void WangLandau::save_state(std::string filename)
{
    Output op;
    op.open(filename.c_str(),"w");
    op<<"# f = "<<gfact<<"\n";
    op<<"# emin = "<<emin<<"\n";
    op<<"# emax = "<<emax<<"\n";
    op<<"# n_bins = "<<nbins<<"\n";
    op<<"# ebin = "<<ebin<<"\n";
    op<<"# flatness threshold = "<<converr<<"\n";
    for (size_t i=0;i<gvals.size();++i) {
        op<<(emin+ebin*(0.5+i))<<"\t"<<his[i]<<"\t"<<gvals[i]<<"\n";
    }
    op.close();
}

void WangLandau::read_g(std::string filename)
{
    if (!prf_utils::STestFile(filename)) return;
    std::ifstream fin(filename.c_str());
    std::string line;
    std::deque<std::string> lines;
    while (getline(fin,line)) {
        line=prf_utils::trim_str(line);
        if (line[0]=='#') {
            std::deque<std::string> tokens;
            prf_utils::split(line,tokens);
            if (tokens.size()<4) continue;
            if (tokens[1]=="emin" && !usr_energy_range) {
                emin=strtod(tokens[3].c_str(),NULL);
            } else if (tokens[1]=="emax" && !usr_energy_range) {
                emax=strtod(tokens[3].c_str(),NULL);
            } else if (tokens[1]=="n_bins" && !usr_n_bins) {
                nbins=atoi(tokens[3].c_str());
            }
        } else {
            lines.push_back(line);
        }
    }
    init_his();
    init_g();
    for (size_t i=0;i<lines.size();++i) {
        line=lines[i];
        std::istringstream ssin(line);
        double e,h,g;
        ssin>>e>>h>>g;
        if (e>=emax) continue;
        if (e<emin) extend_till(e);
        gvals[bin_of(e)]=g;
    }
    fin.close();
    gvalsread=true;
}

int WangLandau::parseCommand(InstructionString s)
{
    if (s.head()=="energy_range") {
        emin=strtod(s.tail().part(0).c_str(),NULL);
        emax=strtod(s.tail().part(1).c_str(),NULL);
        usr_energy_range=true;
    } else if (s.head()=="number_of_bins") {
        n_bins(atoi(s.tail().str().c_str()));
        usr_n_bins=true;
    } else if (s.head()=="flatness_cutoff") {
        flatness_cutoff(strtod(s.tail().str().c_str(),NULL));
    } else if (s.head()=="use_htmc") use_htmc(s.tail().str()=="on");
    else if (s.head()=="initial_wl_factor")
        init_g_factor(strtod(s.tail().str().c_str(),NULL));
    else if (s.head()=="read_g_from") read_g(s.tail().str());
    return 1;
}

int WangLandau::check_bins()
{
    if (nbins!=prevbins) {
        prf::cout<<"There are new bins in the WL histograms. "
                <<"Keeping estimate of g obtained so far in the old bins. \n";
        prevbins=n_bins();
        return 0;
    }
    return 1;
}

int WangLandau::update_g_factor()
{
    if (converged()) {
        char wlstats[50];
        sprintf(wlstats,"%swlstats_%d",output_prefix.c_str(),itmp);
        save_state(wlstats);

        gfact=sqrt(gfact);
        lngf=log(gfact);

        reset_his();
        ++itmp;
    }
    return itmp;
}

void WangLandau::print_setup()
{
    prf::cout<<"\n1 Monte Carlo Cycle = "<<CycleLength()
    <<" Elementary Monte Carlo Steps or Updates\n\n";
    uph.print_setup();
    prf::cout<<"Initial energy range :"<<e_lower_limit()
            <<" to "<<e_upper_limit()<<" with "<<n_bins()
            <<" bins.\n";
    prf::cout<<"The histogram of energy will be regarded flat when "
            <<"the population of the least populated bin is greater "
            <<"than "<<(1-flatness_cutoff())<<" times the population "
            <<"of the most populated bin.\n";
    prf::cout<<"The WL multiplicative factor will be updated as f-->sqrt(f) "
            <<ntmp<<" times.\n";
}

std::string WangLandau::ConfSignature()
{
    std::ostringstream ost;
    ost<<"simulation_method "<<descr<<"\n";
    ost<<"ntmp "<<ntmp<<"\n";
    return ost.str();
}
