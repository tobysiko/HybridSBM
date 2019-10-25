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

#include "PopBase.hh"
#include "../Aux/profasi_io.hh"
#include "../Aux/fileutils.hh"
#include <deque>

using std::string;
using std::list;
using std::deque;

using namespace prf_utils;

namespace prf
{
    SelRes::SelRes()
    {
        mdl=0; chn=0;indx_nat=0;
        indx_str="";nam="";
    }

    SelRes::SelRes(int mdlno,int chno,int i)
    {
        mdl=mdlno;chn=chno;indx_nat=i;
    }

    SelRes::SelRes(const SelRes &pr2)
    {
        mdl=pr2.mdl;chn=pr2.chn;indx_nat=pr2.indx_nat;
        indx_str=pr2.indx_str; nam=pr2.nam;
    }

    SelRes::~SelRes() {}

    SelRes & SelRes::operator=(const SelRes &pr2)
    {
        if (this!=&pr2) {
            mdl=pr2.mdl;chn=pr2.chn;indx_nat=pr2.indx_nat;
            indx_str=pr2.indx_str; nam=pr2.nam;
        }

        return *this;
    }

    bool SelRes::operator==(SelRes pr2)
    {
        return mdl==pr2.mdl && chn==pr2.chn && indx_nat==pr2.indx_nat &&
               indx_str==pr2.indx_str && nam==pr2.nam;
    }

    void SelRes::setval(int mdlno,int chno,int i, string ist, string nm)
    {
        mdl=mdlno;chn=chno;indx_nat=i; indx_str=ist; nam=nm;
    }

    PopBase::PopBase() {nc=0;curmodel=0;nmodels=1;}

    PopBase::~PopBase() {}

    std::string PopBase::chain_name(int i) const {return "?";}

    int PopBase::chain_number(std::string chnm) const
    {
        if (!chnm.empty() and chnm[0]=='#') {
            std::string numstr=chnm.substr(1);
            int i=atoi(numstr.c_str());
            if (i>=0 and i<nc) return i;
            else return 0;
        }
        for (int i=0;i<nc;++i) {
            if (chain_name(i)==chnm) return i;
        }

        return -1;
    }

    int PopBase::num_grp(int ich) const {return 0;}

    int PopBase::index_of_grp(std::string ires,int ich){return -1;}

    string PopBase::grp_name(int ires,int ic) {return "NUL";}

    string PopBase::str_index(int ires,int ic)
    {
        char num[16];
        sprintf(num,"%d",ires);
        return string(num);
    }

    int PopBase::descriptors(std::list<SelRes> &slc,
                             std::list<AtomDescriptor> &des){des.clear();return 0;}

    int PopBase::export_shape(std::vector<int> &vct,Shape &shp){return 0;}

    int PopBase::select(string chid, string s1, string s2,list<SelRes> &slc)
    {
// slc.clear();
        int nselect=0,r1=0,r2=0;
        int chnum=-1;
        if (chid=="*") chnum=nc; else chnum=chain_number(chid);

        if (chnum<0) {
            prf::cerr<<"Selection "<<s1<<" to "<<s2<<": Chain \""
            <<chid<<"\" not found.\n";
            return nselect;
        } else if (chnum==nc) {
            for (int j=0;j<nc;++j) {
                if ((r1=index_of_grp(s1,j))==num_grp(j)) r1=0;

                if ((r2=index_of_grp(s2,j))==num_grp(j)) r2=num_grp(j)-1;

                for (int k=r1;k<=r2;++k) nselect+=add_selection(j,k,slc);
            }
        } else {
            if ((r1=index_of_grp(s1,chnum))==num_grp(chnum)) r1=0;

            if ((r2=index_of_grp(s2,chnum))==num_grp(chnum)) r2=num_grp(chnum)-1;

            for (int k=r1;k<=r2;++k) {
                nselect+=add_selection(chnum,k,slc);
            }
        }

        return nselect;
    }

    int PopBase::add_selection(int ich,int ires,list<SelRes> &slc)
    {
        SelRes s;
        s.mdl=curmodel;
        s.chn=ich;
        s.indx_nat=ires;
        s.indx_str=str_index(ires,ich);
        s.nam=grp_name(ires,ich);
        slc.push_back(s);
        return 1;
    }

    int PopBase::mk_selection(string slcstr, list<SelRes> &slc)
    {
        Logger blog;
        slc.clear();
        int nselect=0;

        if (slcstr.empty()) return nselect;;

        deque<string> chunks,rng;

        split_str<deque<string> >(slcstr,':',chunks);

        int modelno=atoi(chunks.front().c_str());

        set_model(modelno-1);

        blog(40)<<"Trying to make a selection of "<<slcstr<<"\n";

        blog(40)<<"set model number to "<<modelno-1<<"\n";

        chunks.pop_front();

        bool nodetails=true;

        for (deque<string>::iterator i=chunks.begin();i!=chunks.end();++i) {
            if (!(i->empty())) {nodetails=false; break;}
        }

        if (nodetails) nselect+=select("*","undef","undef",slc);
        else {
            for (deque<string>::iterator i=chunks.begin();i!=chunks.end();++i) {
                if (i->empty()) continue;

                blog(40)<<"Evaluating selection "<<(*i)<<"\n";

                split_str<deque<string> >(*i,',',rng);

                string chid=rng.front();

                string res1("undef"),res2("undef");

                if (rng.size()>1) res1=rng[1];

                if (rng.size()>2) res2=rng[2];

                nselect+=select(chid,res1,res2,slc);
            }
        }

        return nselect;
    }

    void PopBase::set_model(int i)
    {
        if (i<nmodels&&i>=0) curmodel=i;
    }
}

