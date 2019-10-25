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

#include "UpdatesHandler.hh"
#include <fstream>
#include <sstream>
#include "../Aux/fileutils.hh"
#include "../Aux/prf_xml.hh"

using namespace prf_utils;

bool UpdatesHandler::is_update_name(std::string nm)
{
    return (nm=="Rot" or nm=="Pivot" or nm=="BGS"
            or nm=="Rotation" or nm=="Translation");
}

Update * UpdatesHandler::new_update(std::string upname)
{
    Update *ans=NULL;
    if (upname=="Rot") ans=new Rot();
    else if (upname=="Pivot") ans=new Pivot();
    else if (upname=="BGS") ans=new BGS();
    else if (upname=="Translation") ans=new Translation();
    else if (upname=="Rotation") ans=new Rotation();
    return ans;
}

Update * UpdatesHandler::used_update(std::string upnm)
{
    Update *ans=NULL;
    for (size_t i=0;i<used.size();++i) {
        if (used[i]->Name()==upnm) {
            ans=used[i];
            break;
        }
    }
    return ans;
}

bool UpdatesHandler::have_update(std::string upname)
{
    bool ans=false;
    for (size_t i=0;i<used.size();++i) {
        if (used[i]->Name()==upname) {
            ans=true;
            break;
        }
    }
    return ans;
}

int UpdatesHandler::use_update(std::string upname)
{
    int added=0;
    if (!have_update(upname)) {
        Update *updt=new_update(upname);
        if (updt!=NULL) {
            used.push_back(updt);
            updt->set_RandomNumberGenerator(rng);
            added=1;
        }
    }
    return added;
}

int UpdatesHandler::use_update(Update *updt)
{
    int added=0;
    if (updt!=NULL and (!have_update(updt->Name()))) {
        used.push_back(updt);
        added=1;
    }
    return added;
}

int UpdatesHandler::skip_update(std::string upname)
{
    int skipped=0;
    if (have_update(upname)) {
        std::vector<Update *> bkp=used;
        used.clear();
        for (size_t i=0;i<bkp.size();++i) {
            if (bkp[i]->Name()!=upname)
                used.push_back(bkp[i]);
            else delete (bkp[i]);
        }
        skipped=1;
    }
    return 1;
}

void UpdatesHandler::init()
{
    for (size_t i=0;i<used.size();++i) used[i]->init();
}

/**
  \page settings_mc Commands for the Monte Carlo
  \section settings_updates Conformational updates
  \li \b --update_prob_file or \b -upfile : Name of the file with probabilities
  for various updates. The format for update probability specification is
  described in \ref updateprob.
  \li \b --config_update or \b -cup : Configure updates using information in
  an XML file. See \ref config_update for details.

  There are also switches for each conformational update, Pivot, BGS, Rot,
  Rotation and Translation. To turn off Pivot for a simulation, use
  for example, --no-Pivot. To force the use of rigid body updates in
  single chain simulations, you can use the switches --Rotation --Translation.
  To specify a switch in the settings file, you write, for instance,
  "Rotation on".
  */

UpdatesHandler::UpdatesHandler() : HandlerBase()
{
    updtprobfile="updateprob.dat";
    ntmp=1;
    uplog=10;
    par.option("update_prob_file", "upfile",1,"(File with update probabilities)");
    par.option("config_updates","cup",1,
               "(configure updates using info in an XML file)");
    par.new_switch("Rot","Rot",true,"(default)");
    par.new_switch("Pivot","Pivot",true,"(default)");
    par.new_switch("BGS","BGS",true,"(default)");
    par.new_switch("Rotation","Rotation",true,"(default, for multi-chain)");
    par.new_switch("Translation","Translation",true,"(default for multi-chain)");
}

UpdatesHandler::~UpdatesHandler()
{
    for (size_t i=0;i<used.size();++i) {
        if (used[i]) {
            delete (used[i]);
            used[i]=NULL;
        }
    }
}

void UpdatesHandler::set_n_temps(size_t nt)
{
    ntmp=nt;
    num_calls.allocate(used.size(), ntmp);
    num_acc_calls.allocate(used.size(), ntmp);
    probs.allocate(used.size(),ntmp);
}

void UpdatesHandler::autoSelect(Population *popl)
{
    p=popl;
    if (p->NumberOfChains()==0) {
        prf::cerr<<"UpdatesHandler> Auto select updates called \n"
                <<"with population of zero size. This may be harmless, \n"
                <<"but it is unlikely that it was intended. It is recommended\n"
                <<"that you check the cause...\n";
        return;
    }

    use_update("BGS");
    use_update("Rot");
    use_update("Pivot");

    if (p->NumberOfChains()>1) {
        use_update("Rotation");
        use_update("Translation");
    }
    set_population(p);
}

void UpdatesHandler::set_population(Population * popl)
{
    p=popl;
    for (size_t i=0;i<used.size();++i) used[i]->connect(p);
}

int UpdatesHandler::parseCommand(InstructionString s)
{
    if (is_update_name(s.head())) {
        if (s.tail().str()=="on") use_update(s.head());
        else if (s.tail().str()=="off") skip_update(s.head());
        else if (s.tail().head()=="set") {
            Update *up=used_update(s.head());
            if (up!=NULL) up->set(s.tail().tail());
        }
    } else if (s.head()=="update_prob_file") updtprobfile=s.tail().str();
    else if (s.head()=="config_updates" or s.head()=="configure_updates") {
        std::string filename=s.part(1);
        prf_xml::XML_Node *upc=prf_xml::get_xml_tree(filename);
        if (upc!=NULL) {
            upc->interpret_formatted_data();
            for (size_t c=0;c<upc->n_children();++c) {
                prf_xml::XML_Node *nd=upc->child(c);
                if (nd->name()=="all_updates") {
                    for (size_t i=0;i<used.size();++i) used[i]->configure(nd);
                } else {
                    Update *u=used_update(nd->name());
                    if (u!=NULL) {
                        u->configure(nd);
                    } else {
                        prf::cerr<<"No update of name "<<nd->name()<<"\n";
                    }
                }
            }
        }
    }

    return 1;
}

void UpdatesHandler::assign_probs()
{
    if (ReadFile(updtprobfile)==0) auto_assign_probs();
}

void UpdatesHandler::auto_assign_probs()
{
    Logger blog;

    if (ntmp<=0) {
        prf::cerr<<"Error! Asked to automatically set up probabilities"
        <<" for "<<ntmp<<" temperatures! Program aborts.\n";
        exit(1);
    }

    int nup= used.size();

    blog(10)<<"Using automatic probability setup for "<<ntmp
    <<" temperatures and "<<nup<<" updates.\n";
    probs.allocate(nup,ntmp);
    //setting probabilities with information about degrees of freedom
    double nc=p->NumberOfChains(),nrtdof=0,nbbdof=0,
                                            nrgdof=0, ntotdof=0;

    for (int i=0;i<nc;++i) {
        nrtdof+=p->Chain(i)->numRTdof();
        nbbdof+=p->Chain(i)->numBBdof();
        nrgdof+=6;
    }

    ntotdof=nrtdof+nbbdof+nrgdof;

    double fracrot,fracbb,fracrg;
    fracrot=nrtdof/ntotdof;
    fracbb=nbbdof/ntotdof;
    fracrg=nrgdof/ntotdof;
    blog(uplog)<<"Fractions of rot, bb, and rigid dof are "<<fracrot
    <<", "<<fracbb<<", "<<fracrg<<"\n";
    int ncu=0;

    for (int i=0;i<nup;++i) {
        if (used[i]->rigid_chain_update()) ++ncu;
    }

    for (int i=0;i<nup;++i) {
        if (used[i]->rigid_chain_update()) {
            for (size_t j=0;j<ntmp;++j) probs[i][j]=fracrg/ncu;
        }
    }

    ncu=0;

    for (int i=0;i<nup;++i) {
        if (used[i]->sidechain_update()) ++ncu;
    }

    for (int i=0;i<nup;++i) {
        if (used[i]->sidechain_update()) {
            for (size_t j=0;j<ntmp;++j) probs[i][j]=fracrot/ncu;
        }
    }

    ncu=0;

    int locu=0;

    for (int i=0;i<nup;++i) {
        if (used[i]->backbone_update()) {
            ++ncu;

            if (used[i]->local_update()) ++locu;
        }
    }

    if (ntmp==1) {
        for (int i=0;i<nup;++i) {
            if (used[i]->backbone_update()) {
                probs[i][0]=fracbb/ncu;
            }
        }
    } else {
        for (int i=0;i<nup;++i) {
            if (used[i]->backbone_update()) {
                for (size_t j=0;j<ntmp;++j) {
                    double fraclocal=0;

                    if (j>=(ntmp/2) && locu!=0)
                        fraclocal=0.2+0.6*(1+j-(ntmp/2))/(ntmp-ntmp/2);

                    if (used[i]->local_update()) {
                        probs[i][j]=fracbb*fraclocal/locu;
                    } else {
                        probs[i][j]=fracbb*(1-fraclocal)/(ncu-locu);
                    }
                }
            }
        }
    }
    normalize();

    blog(10)<<"Based on the information about the population and the \n"
    <<"updates used, the following probability matrix was calculated\n"
    <<"for update probabilities...\n";
    blog<<"probability values ...\n";

    for (int i=0;i<nup;++i) {
        blog<<(used[i]->Name()).c_str()<<": ";

        for (size_t j=0;j<ntmp;++j) blog <<probs[i][j]<<"\t";

        blog<<"\n";
    }

    blog<<"\n";
}

int UpdatesHandler::ReadFile(std::string fl)
{
    Logger blog;
    std::vector<std::string> upnames,vs;
    std::vector<int> torder,uorder;
    bool col_t_row_u=true;

    if (!STestFile_r(fl.c_str())) {
        blog(3)<<"File "<<fl<<" with update probabilities could not be opened. "
                <<"Update probabilities will not be read in from a file.\n";
        return 0;
    }
    std::ifstream fin(fl.c_str());
    std::list<std::string> lines;
    std::string line;

    while (getline(fin,line)) if (!line.empty()) lines.push_back(line);

    fin.close();

    int ilin=0,l=0;

    int nupdates=used.size();

    probs.allocate(nupdates,ntmp);

    std::list<std::vector<double> > plines;

    for (std::list<std::string>::iterator it=lines.begin();it!=lines.end();++it) {
        blog(uplog)<<"line "<<ilin++<<": "<<*it<<"\n";
        vs.clear();
        split<std::vector<std::string> >(*it,vs);

        if (vs.front()=="#Updates") {
            uorder.resize(vs.size()-1,-1);

            for (size_t j=1;j<vs.size();++j) {
                for (int k=0;k<nupdates;++k) {
                    if (used[k]->Name()==vs[j]) {
                        uorder[j-1]=k;
                        upnames.push_back(vs[j]);
                    }
                }

                if (uorder[j-1]==-1) {
                    blog<<"Update "<<vs[j]<<" is not being used.\n";
                    blog<<"Data corresponding to this update from the file "
                    <<"will be ignored.\n";
                }
            }
        } else if (vs.front()=="#T_indices") {
            std::list<int> tinl;

            for (size_t j=1;j<vs.size();++j) {
                if (vs[j]!="-") tinl.push_back(atoi(vs[j].c_str()));
                else {
                    int lowlim=0,upplim=0;

                    if (j>1) lowlim=atoi(vs[j-1].c_str());

                    if (j<vs.size()-1) upplim=atoi(vs[j+1].c_str());

                    if (lowlim<=upplim) {
                        blog<<"Trying to fill from "<<lowlim<<" to "<<upplim<<"\n";

                        for (int k=lowlim+1;k<upplim;++k) tinl.push_back(k);
                    } else {
                        blog<<"Trying to fill from "<<lowlim<<" to "<<upplim<<"\n";

                        for (int k=lowlim-1;k>upplim;--k) tinl.push_back(k);
                    }
                }
            }

            torder.resize(tinl.size(),0);

            std::list<int>::iterator nlt=tinl.begin();

            for (size_t a=0;a<torder.size();++a,++nlt) torder[a]=*nlt;
        } else if (vs[0]=="#Columns") {
            if (vs[1]=="Updates") col_t_row_u=false;
        } else {
            std::vector<double> pline(vs.size(),0.0);

            for (size_t k=0;k<vs.size();++k)
                pline[k]=strtod(vs[k].c_str(),NULL);

            plines.push_back(pline);
        }
    }

    blog<<"Determining maximum number of columns\n";

    size_t maxcols=0;

    for (std::list<std::vector<double> >::iterator it=plines.begin();
         it!=plines.end();++it) {
        maxcols=std::max(maxcols, it->size());
    }

    blog<<"Smoothing out, so that all rows have same number of columns\n";

    for (std::list<std::vector<double> >::iterator it=plines.begin();
         it!=plines.end();++it) {
        if (it->size()<(size_t) maxcols) {
            for (size_t k=it->size();k<maxcols;++k) it->push_back(it->back());
        }
    }

    l=0;

    Matrix<double> fprob(plines.size(),maxcols);

    for (std::list<std::vector<double> >::iterator it=plines.begin();
         it!=plines.end();++it) {
        for (size_t i=0;i<it->size();++i) fprob[l][i]=(*it)[i];

        ++l;
    }

    if (col_t_row_u) {
        for (size_t i=0;i<uorder.size();++i) {
            if (uorder[i]==-1) continue;

            for (size_t j=0;j<std::min(maxcols,ntmp);++j) {
                probs[uorder[i]][torder[j]]=fprob[i][j] ;
            }

            if (maxcols<ntmp) {
                for (size_t j=maxcols;j<ntmp;++j)
                    probs[uorder[i]][torder[j]]=fprob[i][maxcols-1];
            }
        }
    } else {
        for (size_t i=0;i<uorder.size();++i) {
            if (uorder[i]==-1) continue;

            for (size_t j=0;j<std::min(fprob.dim1(),ntmp);++j) {
                probs[uorder[i]][torder[j]]=fprob[j][i] ;
            }

            if (fprob.dim1()<ntmp) {
                for (size_t j=fprob.dim1();j<ntmp;++j)
                    probs[uorder[i]][torder[j]]=fprob[fprob.dim1()-1][i];
            }
        }
    }

    blog<<"This file contains probabilities for the updates...\n";

    for (size_t i=0;i<upnames.size();++i) blog<<upnames[i]<<"\n";

    blog<<"The order of temperature indices in the file is ...\n";

    for (size_t i=0;i<torder.size();++i) blog<<torder[i]<<"\n";

    blog<<"The different columns represent different ";

    if (col_t_row_u) blog<<"temperatures\n";
    else blog<<"updates\n";

    normalize();

    blog(10)<<"In increasing order of temperature indices, the read "
    <<"probability values (after normalization) are \n\n";

    for (int iup=0;iup<nupdates;++iup) {
        blog<<used[iup]->Name()<<": \n";

        for (size_t jtmp=0;jtmp<ntmp;++jtmp) {
            blog<<probs[iup][jtmp]<<"  ";
        }

        blog<<"\n";
    }
    return 1;
}

void UpdatesHandler::normalize()
{
    double nrom=0;

    for (size_t j=0;j<probs.dim2();++j) {
        nrom=0;

        for (size_t i=0;i<probs.dim1();++i) nrom+=probs[i][j];

        for (size_t i=0;i<probs.dim1();++i) probs[i][j]/=nrom;
    }
}

void UpdatesHandler::print_updateprobs(Output &op)
{
    op<<"****** Updates and their probabilities ******\n";
    if (ntmp>1) {
        op<<"Different rows correspond to different temperature "
            <<"indices\n";
        op<<"Index\t";
    }

    for (size_t i=0;i<used.size();++i) {
        op<<used[i]->Name()<<"\t\t";
    }
    op<<"\n";
    for (size_t j=0;j<probs.dim2();++j) {
        if (ntmp>1) op<<j<<"\t";
        for (size_t i=0;i<used.size();++i) {
            op<<probs[i][j]<<"\t";
        }
        op<<"\n";
    }
}

Update * UpdatesHandler::perform_update(int itmp)
{
    double r=rng->shoot();
    int i=0;

    while ((r-=probs(i,itmp))>0) ++i;
    active_update=i;
    used[i]->perform(itmp);
    tindex=itmp;
    num_calls[i][itmp]++;
    return used[i];
}

void UpdatesHandler::accept_update()
{
    used[active_update]->accept();
    num_acc_calls[active_update][tindex]++;
}

void UpdatesHandler::reject_update()
{
    used[active_update]->revert();
}

void UpdatesHandler::output_statistics(std::string filename)
{
    Output fp(filename.c_str());
    fp<<"Acceptance information for various updates...\n";
    for (size_t i=0;i<used.size();++i) {
        double ncalls=0,nacc=0;
        fp<<"############### "<<used[i]->Name()<<" ###############\n";
        fp<<"Temperature index  |  Total calls | Accepted calls | Acceptance fraction\n";
        for (size_t j=0;j<ntmp;++j) {
            ncalls+=num_calls[i][j];
            nacc+=num_acc_calls[i][j];
            fp<<j<<"\t|\t"<<(num_calls[i][j])<<"\t|\t"
                    <<num_acc_calls[i][j]<<"\t|\t"
                    <<((num_calls[i][j]!=0)?((double) num_acc_calls[i][j])/(num_calls[i][j]):0)
                    <<"\n";
        }
        fp<<"Over all, "<<ncalls<<" calls, "<<nacc<<" accepted. \n\n";
    }

    fp.close();
}

void UpdatesHandler::reset_statistics()
{
    num_calls*=0;
    num_acc_calls*=0;
}

void UpdatesHandler::setBeta(double bt)
{
}

void UpdatesHandler::print_setup()
{
    std::string st;
    for (size_t i=0;i<used.size();++i) used[i]->print_setup(st);
    if (!st.empty()) {
        highlight("Modifications of update behaviour");
        prf::cout<<st<<"\n";
        highlight();
    }
}
