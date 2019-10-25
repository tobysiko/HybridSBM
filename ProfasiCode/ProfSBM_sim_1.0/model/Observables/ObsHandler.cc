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

#include "ObsHandler.hh"
#include "../Aux/fileutils.hh"
#include "Rg.hh"
#include "AtomDistance.hh"
#include "ProteinRMSD.hh"
#include "NativenessQ.hh"
#include "NativenessQnoAln.hh"
#include "RCAngleRMSD.hh"
#include "ContactOrder.hh"
#include "SecStruct.hh"
#include "HelixSegment.hh"
#include "ContactMap.hh"
#include "OligoOrient.hh"
#include "../Energy/FFHandler.hh"
#include "PeriodicObs.hh"
#include <sstream>
#include <cmath>

using namespace prf;

using namespace prf_utils;
using std::vector;
using std::string;
using std::map;
using std::deque;
using std::find;
using std::ifstream;

ObsHandler::ObsHandler()
{
    myprefx="tmp_";
    avgfile="averages";
    rtfile="rt";
    dirc=".";
    nobs=curst=curtime=0;
    nstates=1;
    obs.clear();
    obsno.clear();
    hs2d.clear();
    hs2dreq1.clear();
    hs2dreq2.clear();
    //Set up program switches
    swtch[string("histograms")]=true;
    swtch[string("hist_resume")]=false;
    //None of the switches set here were specified by the user...

    for (map<string,bool>::iterator it=swtch.begin();it!=swtch.end();++it) {
        usrspc[it->first]=false;
    }
    par.option("new_obs","new_obs",3,"(Observable type, alias and options inside enclosing quotes)");
    par.option("set_obs","set_obs",2,"(Observable alias and options inside enclosing quotes)");
    par.option("make_2d_his","make_2d_his",2,"(Aliases of the two observables)");
    par.new_switch("histograms","mkhis",true,
                    "(make histograms during the run)");
    par.new_switch("hist_resume","hist_resume",false,"(Restore histogram state while resuming runs)");
}

ObsHandler::~ObsHandler()
{
    rt.close();

    for (size_t i=0;i<hs2d.size();++i) {
        if (hs2d[i]) delete hs2d[i];
    }

    for (size_t i=0;i<myobs.size();++i) {
        if (myobs[i]) {
            delete myobs[i];
            myobs[i]=NULL;
        }
    }
}

void ObsHandler::setDir(string write_dir)
{
    dirc=write_dir;
    setPrefix(write_dir+string("/"));
}

void ObsHandler::setPrefix(string prfx)
{
    myprefx=prfx;
    avgfile=myprefx+avgfile;
    rtfile=myprefx+rtfile;

    for (size_t i=0;i<obs.size();++i) obs[i]->output_prefix(myprefx);
}

Observable * ObsHandler::get_obs(std::string nm)
{
    Observable *ans=NULL;

    for (size_t i=0;i<obs.size();++i) {
        if (obs[i]->Name()==nm) {ans=obs[i];break;}
    }

    return ans;
}

void ObsHandler::sample(int current_time,int current_state)
{
    curtime=current_time;
    curst=current_state;

    for (size_t i=0;i<obs.size();++i) {
        obs[i]->refresh();
        obs[i]->avg_fill(curst);
        obs[i]->his_fill(curst);
    }

    for (size_t i=0;i<hs2dobs1.size();++i)
        hs2d[nstates*i+curst]->put(hs2dobs1[i]->Value(),hs2dobs2[i]->Value());
}

void ObsHandler::enableStatistics()
{
    for (size_t i=0;i<obs.size();++i) {
        obs[i]->enable_stats();
    }
}

void ObsHandler::disableStatistics()
{
    for (size_t i=0;i<obs.size();++i) {
        obs[i]->disable_stats();
    }
}

void ObsHandler::writeAverages()
{
    Output avg(avgfile.c_str(),"w");
    avg<<"# Mean and standard deviation of various observables \n";
    avg<<"# "<<stdesg<<"-index\tmean\tstandard-deviation  \n";

    for (size_t i=0;i<obs.size();++i) obs[i]->avg_write(avg);

    avg.close();
}

void ObsHandler::writeHistograms()
{
    for (size_t i=0;i<obs.size();++i) {
        obs[i]->his_save();
    }

    save2DHis();
}

void ObsHandler::writeRTSnapshot()
{
    rt<<curtime<<"  "<<curst;

    for (size_t i=0;i<obs.size();++i) {
    	//std::cout<<"obs name:"<<obs[i]->Name()<<"\n";
    	obs[i]->write_snapshot(rt);
    }

    rt<<"\n";

    rt.flush();
}

void ObsHandler::resetHistograms()
{
    for (size_t i=0;i<obs.size();++i) obs[i]->his_reset();

    for (size_t i=0;i<hs2d.size();++i) hs2d[i]->reset();
}

void ObsHandler::resetAvgData()
{
    for (size_t i=0;i<obs.size();++i) obs[i]->avg_reset();
}

void ObsHandler::track(Observable *ob,string options)
{
    if (ob==NULL) return;
    add_obs_and_opts(ob,ob->Name(),options);
}

void ObsHandler::writeRTKey()
{
    Output rtkey((myprefx+"rtkey").c_str(),"w");
    rtkey<<"Different columns of the rt file have the following information\n"
    <<"Time in MC cycles \n"<<stdesg<<"\n";

    for (size_t i=0;i<obs.size();++i) {
        obs[i]->write_rtkey(rtkey);
    }

    rtkey.close();

    rt.open(rtfile.c_str(),"a");
}

int ObsHandler::set_block_props(int num_states, string state_designator)
{
    nstates=num_states;
    stdesg=state_designator;
    for (size_t i=0;i<obs.size();++i) {
        obs[i]->set_n_temp(num_states);
    }
    return 1;
}

void ObsHandler::fix_dependencies(ObsEnergy * depterm)
{
    Observable * tmpo=get_obs(depterm->obsName());
    if (tmpo!=NULL) depterm->connectObs(tmpo);
    else {
        prf::cerr<<"ObsHandler> Could not satisfy dependency "
                <<depterm->obsName()<<" for dependent object "
                <<depterm->Name()<<".\n";
    }
}

int ObsHandler::initialize()
{
    for (size_t i=0;i<obs.size();++i) {

        if (swtch[string("hist_resume")]) obs[i]->set_his_resume();

        if (swtch[string("histograms")]) obs[i]->make_his(true);

        obs[i]->his_setup();
    }

    for (size_t i=0;i<hs2dreq1.size();++i) {
        make2DHis(hs2dreq1[i],hs2dreq2[i]);
    }

    setup2DHis();

    return 1;
}

void ObsHandler::printConfig()
{
    highlight("Observables");
    prf::cout<<"Run-time history file \"rt\" will contain ...\n\n";
    prf::cout<<"0. Time in MC cycles. \n";
    prf::cout<<"1. Temperature index. \n";

    for (size_t i=0; i<obs.size();++i) {
        prf::cout<< (i+2) <<". ";
        obs[i]->write_rtkey(prf::cout);
    }

    prf::cout<<"\nThe averages file will contain averages for all the "

    <<"Observables in the rt file, as mentioned above.\n\n";

    if (hs2dobs1.size() >0) {
        prf::cout<<"The following 2d (correlation) histograms will be made \n";

        for (size_t i=0;i<hs2dobs1.size();++i) {
            prf::cout<<hs2dobs1[i]->Name() <<" vs ";
            prf::cout<<hs2dobs2[i]->Name() <<" : ("<<nstates<<" histograms)\n";
            prf::cout<<"X-range : ["<<hs2d[nstates*i]->Xmin() <<", "
            <<hs2d[nstates*i]->Xmax() <<"), with "
            <<hs2d[nstates*i]->NXbins() <<" bins.\n";
            prf::cout<<"Y-range : ["<<hs2d[nstates*i]->Ymin() <<", "
            <<hs2d[nstates*i]->Ymax() <<"), with "
            <<hs2d[nstates*i]->NYbins() <<" bins.\n";
        }
    }
}

void ObsHandler::save2DHis(int i, int j)
{
    char filename[100]="";
    sprintf(filename,"%s%s_%s_%d",myprefx.c_str(),hs2dobs1[i]->Name().c_str(),
            hs2dobs2[i]->Name().c_str(),j);
    hs2d[nstates*i+j]->save_state(filename);
}

void ObsHandler::save2DHis()
{
    for (size_t i=0;i<hs2dobs1.size();++i) {
        for (int j=0;j<nstates;++j) {
            save2DHis(i,j);
        }
    }
}

void ObsHandler::setup2DHis()
{
    Logger blog;

    for (size_t i=0;i<hs2dobs1.size();++i) {
        blog(5) <<"adding "<<nstates<<" 2d histograms for "
        <<hs2dobs1[i]->Name() <<" vs "<<hs2dobs2[i]->Name() <<"\n";

        for (int j=0;j<nstates;++j) {
            His2D *newhis = new His2D();
            hs2d.push_back(newhis);
        }

        double xmn,xmx,ymn,ymx;

        int nxbins=100,nybins=100;
        hs2dobs1[i]->rangeEstimate(xmn,xmx);
        hs2dobs2[i]->rangeEstimate(ymn,ymx);
        size_t histindx=0;

        for (;histindx<obs.size();++histindx) {
            if (obs[histindx]==hs2dobs1[i]) break;
        }

        if (histindx!=obs.size()) {
            xmn=obs[histindx]->histogram()->Xmin();
            xmx=obs[histindx]->histogram()->Xmax();
            nxbins= obs[histindx]->histogram()->Nbins();

            if (nxbins>100) nxbins=100;
        }

        histindx=0;

        for (;histindx<obs.size();++histindx) {
            if (obs[histindx]==hs2dobs2[i]) break;
        }

        if (histindx!=obs.size()) {
            ymn=obs[histindx]->histogram()->Xmin();
            ymx=obs[histindx]->histogram()->Xmax();
            nybins= obs[histindx]->histogram()->Nbins();

            if (nybins>100) nybins=100;
        }

        for (int j=0;j<nstates;++j) {
            hs2d[nstates*i+j]->Range(xmn,xmx,ymn,ymx);
            hs2d[nstates*i+j]->NXbins(nxbins);
            hs2d[nstates*i+j]->NYbins(nybins);
            hs2d[nstates*i+j]->init();
        }
    }
}

int ObsHandler::make2DHis(string obsname1, string obsname2)
{
    bool flag1=false,flag2=false;
    int i1=-1, i2=-1;

    for (size_t i=0;i<obs.size();++i) {
        if (obs[i]->Name() ==obsname1) {flag1=true;i1=i;}

        if (obs[i]->Name() ==obsname2) {flag2=true;i2=i;}

        if (flag1 && flag2) break;
    }

    if (!flag1) prf::cerr<<"ObsHandler: "<<obsname1
        <<" is not among the active observables\n";

    if (!flag2) prf::cerr<<"ObsHandler: "<<obsname2
        <<" is not among the active observables\n";

    if (!(flag1 && flag2)) {
        prf::cerr<<"Unable to make a 2D histogram of "
        <<obsname1<<" and "<<obsname2<<"\n";
        return 0;
    }

    if (obsname1 == obsname2) {
        prf::cerr<<"Making a 2D histogram of observable "<<obsname1<<"\n"
        <<"with itself is ill-advised. It contains no useful \n"
        <<"information, and would just slow down the program, \n"
        <<"while eating memory. Refusing to add 2D histogram.\n";
        return 0;
    }

    hs2dobs1.push_back(obs[i1]);

    hs2dobs2.push_back(obs[i2]);

    return 1;
}

int ObsHandler::add_2d_his(InstructionString s)
{
    if (s.n_parts()>=2) {
        hs2dreq1.push_back(s.part(0));
        hs2dreq2.push_back(s.part(1));
    }
    return 1;
}

int ObsHandler::parseCommand(InstructionString s)
{
    int handled=0;
    if (s.head()=="new_obs") {
        make_obs(s.tail());
        handled=1;
    } else if (s.head()=="set_obs") {
        set_obs_props(s.tail());
        handled=1;
    } else if (s.head()=="add_2d_his") {
        add_2d_his(s.tail());
        handled=1;
    } else {
        for (map<string,bool>::iterator it=swtch.begin();it!=swtch.end();++it) {
            if (s.head()==it->first) {
                handled=1;
                if (s.tail().str()=="on") {
                    swtch[s.head()]=true;
                    usrspc[s.head()]=true;
                    break;
                } else {
                    swtch[s.head()]=false;
                    usrspc[s.head()]=true;
                    break;
                }
            }
        }
    }
    return handled;
}

void ObsHandler::make_obs(InstructionString spec)
{
    if (spec.n_parts()==0) return;

    string otype=spec.head(),pars="",name="";

    if (spec.n_parts()>1) name=spec.part(1);

    for (size_t i=2;i<spec.n_parts();++i) pars+=(spec.part(i)+string(" "));

    if (!name.empty()) {
        bool nmavlbl=true;
        for (size_t i=0;nmavlbl&&i<obs.size();++i) if (obs[i]->Name()==name) nmavlbl=false;
        if (!nmavlbl) {
            prf::cerr<<"ObsHandler> Tried to add observable with name "<<name
                    <<". This name is being used by another observable. "
                    <<"Configuring properties of the existing observable "
                    <<"instead of creating a new one.\n";
            return set_obs_props(InstructionString(name+" "+pars));
        }
    }

    Observable *newobs=NULL;

    if (otype=="Rg") newobs=new Rg();
    else if (otype=="AtomDistance") newobs=new AtomDistance();
    else if (otype=="ProteinRMSD") newobs=new ProteinRMSD();
    else if (otype=="NativenessQ") newobs=new NativenessQ();
    else if (otype=="NativenessQnoAln") newobs=new NativenessQnoAln();
    else if (otype=="RCAngleRMSD") newobs=new RCAngleRMSD();
    else if (otype=="ContactOrder") newobs=new ContactOrder();
    else if (otype=="RCBin") newobs=new RCBin();
    else if (otype=="HelixSegment") newobs=new HelixSegment();
    else if (otype=="ContactMap") newobs=new ContactMap();
    else if (otype=="OligoOrient") newobs=new OligoOrient();
    else if (otype=="PeriodicObs") newobs=new PeriodicObs(this);
    else {
        FFHandler tmpffh;
        newobs=tmpffh.new_energy_term(otype);
    }

    if (newobs!=NULL) {
        add_obs_and_opts(newobs,name,pars);
        myobs.push_back(newobs);
    }

}

void ObsHandler::add_obs_and_opts(Observable *tmp, string nm, string pars)
{
    if (tmp==NULL) return;
    tmp->Name(nm);
    tmp->setPopulation(popl);
    tmp->set_n_temp(nstates);
    tmp->output_prefix(myprefx);
    if (!pars.empty()) tmp->set(pars);
    use_obs(tmp,tmp->init_obs());
}

void ObsHandler::set_obs_props(InstructionString spec)
{
    if (spec.n_parts()==0) return;

    Observable *o=get_obs(spec.head());

    if (o==NULL) {
        for (size_t i=0;i<unused.size();++i) {
            if (unused[i]->Name()==spec.head()){
                o=unused[i];
                break;
            }
        }
    }

    if (o!=NULL) {
        o->set(spec.tail().str());
        use_obs(o,o->init_obs());
    }

}

void ObsHandler::use_obs(Observable *o, bool flg)
{
    Logger blog;
    if (flg) {
        if (std::find(obs.begin(),obs.end(),o)==obs.end()) {
            obsno[o->Name()]=obs.size();
            obs.push_back(o);
            blog(10)<<"ObsHandler> activating observable "<<o->Name()<<"\n";
            nobs=obs.size();
        }
        vector<Observable *> tmpunused;
        for (size_t i=0;i<unused.size();++i) {
            if (unused[i]!=o) tmpunused.push_back(unused[i]);
        }
        unused=tmpunused;
    } else {
        if (std::find(unused.begin(),unused.end(),o)==unused.end()) {
            unused.push_back(o);
            blog(10)<<"ObsHandler> deactivating observable "<<o->Name()<<"\n";
        }
        vector<Observable *> tmpused;
        obsno.clear();
        nobs=0;
        for (size_t i=0;i<obs.size();++i) {
            if (obs[i]!=o) {
                tmpused.push_back(obs[i]);
                obsno[obs[i]->Name()]=nobs++;
            }
        }
        obs=tmpused;
    }
}

void ObsHandler::re_init_obs()
{
    for (size_t i=0;i<obs.size();++i) use_obs(obs[i],obs[i]->init_obs());
    for (size_t i=0;i<unused.size();++i) use_obs(unused[i],unused[i]->init_obs());
}

/**
* \page settings_obs Commands for setting up measurements/Observables
If the simulation is to make a certain measurement during the run, there
are two things one has to do: 1. One has to tell the simulation program
what that measurement is. 2. One might want to (have to) fine-tune the
measurement with some additional information.
\li \b --new_obs or \b -new_obs :
Usage example: <i>--new_obs Rg rg_chainA "of_chain 0" </i> \n
The instruction new_obs is used to add a new measurement to the simulation.
The example here demonstrates how to add the radius of gyration of the
first chain in the population as a measurement. "Rg" here is the name
of the observable type. "rg_chainA" is an alias created by the user. This
means, histograms created for this observable will use this alias. The
third argument passed to new_obs is a set of options for the observable.
In this case, we instruct the Rg object to limit the calculations to
chain 0. These options to modulate the behaviour of the measurements
must be passed enclosed in quote marks on the command line. Even if
there is no modulating option, a blank string needs to be passed. In
the settings file, the quote marks as well as the blank string can be
omitted.  Below, you will find a list of various possible measurement
types. Follow the links to see what modulating options each measurement
supports.
\li \b --set_obs or \b -set_obs :
Usage example: <i>--set_obs rg_chainA "his_range 0 15 fixed"</i> \n
The command set_obs allows you to set additional options for an existing
measurement. The first argument here is the same alias used while
adding the observable. The second argument is a set of modulating options.
In this case, we show how to tell the rg_chainA measurement to use a
fixed range histogram between 0 and 15 \AA.
The reason why a separate command exists for this task is : there are
some measurements in ProFASi that the user does not create with "new_obs",
such as the Energy terms. They are always there in the list of measurements.
The user can send additional modulating instructions to such Observables.
Also, it may be cleaner to put some options in a different line in a settings
file.
\li \b --histograms or \b -mkhis : Normally, a histogram is always made for
every measurement in ProFASi. This instruction is most likely to be used
in the negative form, to globally switch off all histograms during the
run. For example: \em --no-histograms. This tells all observables to
make no histograms during the run. Individual observables can still be
asked to make histograms through their own settings through the set_obs
command above, in which case the more specific instruction takes precedence.


When multiple options are passed to an Observable, either with new_obs or
with set_obs, they should be separated by semi-colon (";") marks. There need
not be a semi-colon after the last option. For example, the following
adds a backbone RMSD measurement between residues 20 -- 30 in a PDB file
and the corresponding residues in the simulation:
<i>--new_obs ProteinRMSD bbrmsd "using +BB+CB ; struc1 abc.pdb::A,20,30 ;
struc2: $::A"</i>

\section profasi_obs Available Observables in PROFASI
- \subpage opt_Observable
- \subpage opt_AtomDistance
- \subpage opt_ContactMap
- \subpage opt_ContactOrder
- \subpage opt_NativenessQ
- \subpage opt_OligoOrient
- \subpage opt_ProteinRMSD
- \subpage opt_RCAngleRMSD
- \subpage opt_Rg
- \subpage opt_RCBin
- \subpage opt_PeriodicObs
*/
