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
#include <stdlib.h>
#include "RMSD_Utils.hh"
#include <deque>
#include "AtomRecord.hh"
#include "../Elements/PopBase.hh"

using std::string;
using std::list;
using std::deque;
using std::vector;

using namespace prf_utils;

namespace prf
{
    isAtom::isAtom() {set_id("is_any_atom");}

    isAtom::~isAtom() {}

    bool isAtom::operator ()(AtomDescriptor atd)
    {
        return true;
    }

    isHeavyAtom::isHeavyAtom() {set_id("is_not_H");}

    isHeavyAtom::~isHeavyAtom() {}

    bool isHeavyAtom::operator()(AtomDescriptor atd)
    {
        return atd.atom_type != 'H';
    }

    isBackboneAtom::isBackboneAtom() {set_id("is_on_BB");}

    isBackboneAtom::~isBackboneAtom() {}

    bool isBackboneAtom::operator()(AtomDescriptor atd)
    {
        return atd.atom_label==string(" CA ")
               ||atd.atom_label==string(" C  ")
               ||atd.atom_label==string(" N  ");
    }

    AtomTypeSelector::AtomTypeSelector(char gtyp)
    {
        mytp=gtyp;
        string nn="";
        nn.push_back(gtyp);
        set_id(nn);
    }

    AtomTypeSelector::~AtomTypeSelector() {}

    bool AtomTypeSelector::operator()(AtomDescriptor atd)
    {
        return atd.atom_type == mytp;
    }

    ResSelector::ResSelector(string gtyp)
    {
        mytp=gtyp;
        set_id(gtyp);
    }

    ResSelector::~ResSelector() {}

    bool ResSelector::operator()(AtomDescriptor atd)
    {
        return atd.resnm == mytp;
    }

    LabelSelector::LabelSelector(string lbl)
    {
        label="";
        strct=false;
        int i=-1;

        while ((++i<(int)lbl.size())&&isspace(lbl[i])) ;

        while ((i<(int)lbl.size())&&(!isspace(lbl[i]))) label+=lbl[i++];

        set_id(string("is(")+label+string(")"));
    }

    LabelSelector::~LabelSelector() {}

    bool LabelSelector::operator()(AtomDescriptor atd)
    {
        string des=atd.atom_label;
        string desgn="";
        int i=-1;

        while ((++i<(int)des.size())&&isspace(des[i])) ;

        while ((i<(int)des.size())&&(!isspace(des[i]))) desgn+=des[i++];

        return (desgn==label ||
                (!strct && label.size()>0 &&
                 desgn.find(label.c_str())!=desgn.npos));
    }

    RMSD_Utils::RMSD_Utils()
    {
        log_thres=10;
        clear_filters();
    }

    void RMSD_Utils::clear_filters()
    {
        for (list<SPA *>::iterator i=predlst.begin();i!=predlst.end();++i)
            if ((*i)) delete(*i);

        predlst.clear();

        incl_list.clear();

        excl_list.clear();

        cumul_p.clear();

        cumul_m.clear();
    }

    RMSD_Utils::~RMSD_Utils() {clear_filters();}

    int RMSD_Utils::make_filter(string expr)
    {
        Logger blog(log_thres);
        string compct="";

        for (size_t i=0;i<expr.size();++i) {
            if (!isspace(expr[i])) compct.push_back(expr[i]);
        }

        blog(30*log_thres)<<"make_filter: removed spaces from filter string --> "
                <<compct<<"\n";

        int nsep=0,npsep=0,nmsep=0;
        std::deque<int> sep,incl,excl;
        std::deque<string> incl_str,excl_str;

        for (size_t i=0;i<compct.size();++i) {
            if (compct[i]=='+'||compct[i]=='-') {
                sep.push_back(i);

                if (compct[i]=='+') {incl.push_back(nsep);++npsep;}

                if (compct[i]=='-') {excl.push_back(nsep);++nmsep;}

                ++nsep;
            }
        }

        for (int i=0;i<npsep;++i) {
            int stpos=sep[incl[i]]+1;
            size_t ndpos=incl[i]+1;

            if (ndpos<sep.size()) ndpos=sep[ndpos];
            else ndpos=compct.size();

            if (ndpos-stpos>0)
                incl_str.push_back(compct.substr(stpos,ndpos-stpos));
        }

        for (int i=0;i<nmsep;++i) {
            int stpos=sep[excl[i]]+1;
            size_t ndpos=excl[i]+1;

            if (ndpos<sep.size()) ndpos=sep[ndpos];
            else ndpos=compct.size();

            if (ndpos-stpos>0)
                excl_str.push_back(compct.substr(stpos,ndpos-stpos));
        }

        blog(30*log_thres)<<"make_filter: Number of separators = "<<nsep<<"\n";

        blog(30*log_thres)<<"make_filter: Number of inclusive filters = "
        <<incl_str.size()<<"\n";
        blog(30*log_thres)<<"make_filter: Number of exclusive filters = "
        <<excl_str.size()<<"\n";

        for (size_t i=0;i<incl_str.size();++i) {
            add_filter(incl_str[i],incl_list);
        }

        for (size_t i=0;i<excl_str.size();++i) {
            add_filter(excl_str[i],excl_list);
        }

        if (incl_list.empty()) {
            blog(30*log_thres)<<"RMSD_Utils> Inclusive filter list was empty."
                    <<"Default to is_any_atom\n";
            isAtom *isat=new isAtom();
            predlst.push_back(isat);
            incl_list.push_back(CPA(isat));
        }

        if (excl_list.empty()) {
            blog(30*log_thres)<<"RMSD_Utils> Exclusive filter list is empty. "
                    <<"That's OK.\n";
            SimplePredicate<AtomDescriptor> *sim=
                new SimplePredicate<AtomDescriptor>;
            sim->set_default_state(false);
            predlst.push_back(sim);
            excl_list.push_back(CPA(sim));
        }

        cumul_p.clear();

        list<CPA >::iterator ptt,mtt,jt,it=incl_list.begin();
        jt=it;
        cumul_p.push_back(*it);
        it=cumul_p.begin();

        for (++jt;jt!=incl_list.end();++jt,++it) {
            cumul_p.push_back((*jt)||(*it));
        }

        blog(30*log_thres)<<"RMSD_Utils> Built cumulative inclusion filters.\n";

        ptt=it;
        cumul_m.clear();
        it=excl_list.begin();
        jt=it;
        cumul_m.push_back(*it);
        it=cumul_m.begin();

        for (++jt;jt!=excl_list.end();++jt,++it) {
            cumul_m.push_back((*jt)||(*it));
        }

        blog(30*log_thres)<<"RMSD_Utils> Built cumulative exclusion filters.\n";

        mtt=it;
        mtt->inversion(true);
        total_filter=(*mtt)&&(*ptt);
        blog(log_thres)<<"RMSD_Utils> Created filter : "<<total_filter.id()<<"\n";
        return 1;
    }

    void RMSD_Utils::apply_filter(list<AtomDescriptor> &lst)
    {
        for (list<AtomDescriptor>::iterator i=lst.begin(),j=i;i!=lst.end();) {
            if (!total_filter(*(j=i++))) {
                lst.erase(j);
            }
        }
    }

    void RMSD_Utils::add_filter(string flt,
                                list<CPA > &lst)
    {
        if (flt=="all"||flt=="any"||flt=="everything") {
            isAtom *isat=new isAtom;
            lst.push_back(CPA(isat));
            predlst.push_back(isat);
            return;
        }
        if (flt=="BB"||flt=="backbone"||flt=="Backbone") {
            isBackboneAtom *bbs=new isBackboneAtom;
            lst.push_back(CPA(bbs));
            predlst.push_back(bbs);
            return;
        }

        if (flt=="HV"||flt=="heavyatom"||flt=="Heavyatom") {
            isHeavyAtom *hvs=new isHeavyAtom;
            lst.push_back(CPA(hvs));
            predlst.push_back(hvs);
            return;
        }

        if (flt[0]=='%') {
            ResSelector *rs=new ResSelector(flt.substr(1,flt.size()));
            lst.push_back(CPA(rs));
            predlst.push_back(rs);
            return;
        }

        if (flt[0]=='@') {
            AtomTypeSelector *att=new AtomTypeSelector(flt[1]);
            lst.push_back(CPA(att));
            predlst.push_back(att);
            return;
        }

        LabelSelector *lbs=new LabelSelector(flt);

        lst.push_back(CPA(lbs));
        predlst.push_back(lbs);
    }

    int RMSD_Utils::map_atoms(list<AtomDescriptor> l1,
                              list<AtomDescriptor> l2,
                              vector<int> &v1, vector<int> &v2, string opts)
    {
        Logger blog;
        v1.clear();
        v2.clear();

        if (l1.empty()) {prf::cerr<<"map_atoms: list 1 empty.\n"; return 0;}

        if (l2.empty()) {prf::cerr<<"map_atoms: list 2 empty.\n"; return 0;}

        list<AtomDescriptor>::iterator it1,it2,it3,it4;

        list<AtomDescriptor> match1,match2;
        string rc1,rp1,rc2,rp2;
        it2=it3=it4=l2.begin();
        rc1=rp1=l1.front().ires;
        rc2=rp2=l2.front().ires;

        while (it4!=l2.end()&&it4->ires==rp2) ++it4;

        for (it1=l1.begin();it1!=l1.end();++it1) {
            if ((rc1=it1->ires)!=rp1) {
                it3=it4;
                rp2=rc2;

                while (it4!=l2.end()&&it4->ires==rp2) ++it4;

                if (it4!=l2.end()) rc2=it4->ires;

                rp1=rc1;
            }

            bool matched=false;

            for (it2=it3;it2!=it4;++it2) {
                if (it1->corresponds_to(*it2 , opts)) {
                    if (it2==it3) ++it3;

                    match1.push_back(*it1);

                    match2.splice(match2.end(),l2,it2);

                    matched=true;
                }
            }

            if (!matched) {
                blog(log_thres)<<"RMSD_Utils> No match found for atom "
                        <<it1->short_info()
                        <<" in list2\n";
            } else {
                blog(5*log_thres)<<"RMSD_Utils> Matched "<<it1->short_info()
                        <<"\nto\n"<<it2->short_info()<<"\n";
            }
        }

        for (it2=l2.begin();it2!=l2.end();++it2) {
            blog(log_thres)<<"RMSD_Utils> No match found for atom "<<it2->short_info()
            <<" in list1\n";
        }

        if (match1.size()!=match2.size()) {
            prf::cerr<<"RMSD_Utils> Error in atom matching algorithm. Unequal sized "
            <<"matched-lists created from the given lists.\n";
            exit(1);
        }

        v1.resize(match1.size());

        v2.resize(match1.size());
        int i=0;

        for (it1=match1.begin(),it2=match2.begin();
             it1!=match1.end();++it1,++it2,++i) {
            v1[i]=it1->int_label;
            v2[i]=it2->int_label;
        }

        return 1;
    }

    void RMSD_Utils::alt_coord_policy(string pol)
    {
        mypolicy=pol;
    }

    void RMSD_Utils::alt_coord_parse(list<AtomDescriptor> &lst)
    {
        list<AtomDescriptor>::iterator it,jt,kt,st,qt;
        list<AtomDescriptor> tmp;

        for (it=lst.begin();it!=lst.end();) {
            if (it->ialt==' ') {++it;continue;}

            for (kt=it;kt!=lst.end() && kt->atom_label == it->atom_label
                 && kt->iresrel == it->iresrel; ++kt) ;

            tmp.splice(tmp.end(),lst,it,kt);

            st=tmp.begin();

            if (mypolicy=="occupancy") {
                string mxocc="0.00";

                for (qt=tmp.begin();qt!=tmp.end();++qt) {
                    if (qt->occ>mxocc) {st=qt;mxocc=qt->occ;}
                }
            } else if (mypolicy=="last") {
                st=tmp.end();
                if (st!=tmp.begin()) --st;
            }

            st->ialt=' ';

            lst.splice(kt,tmp,st);
            it=kt;
            tmp.clear();
        }
    }

    int RMSD_Utils::score(vector<string> &v1, vector<string> &v2, size_t i1, size_t i2)
    {
        int ans=0;
        while (i1<v1.size() && i2<v2.size() && v1[i1]==v2[i2]) {
            ++ans;
            ++i1;
            ++i2;
        }
        return ans;
    }
    void RMSD_Utils::align_segs(list<SelRes>::iterator i1, list<SelRes>::iterator i2,
                                list<SelRes>::iterator j1, list<SelRes>::iterator j2,
                                list<SelRes> &resi, list<SelRes> &resj)
    {
        vector<string> sq1,sq2;
        for (list<SelRes>::iterator k=i1;k!=i2;++k) sq1.push_back(k->nam);
        for (list<SelRes>::iterator k=j1;k!=j2;++k) sq2.push_back(k->nam);
        size_t smaller = sq1.size()<sq2.size() ? sq1.size() : sq2.size();
        if (smaller<3) return;

        int mxscr=-1,scr=0;
        size_t maxi=0, maxj=0;
        for (size_t i=0;i<sq1.size();++i) {
            for (size_t j=0;j<sq2.size();++j) {
                scr=score(sq1,sq2,i,j);
                if (scr>mxscr) {
                    maxi=i;
                    maxj=j;
                    mxscr=scr;
                }
            }
        }
        if (mxscr<3) return;
        list<SelRes>::iterator maxit=i1,maxjt=j1;
        for (size_t i=0;i<maxi;++i) ++maxit;
        for (size_t j=0;j<maxj;++j) ++maxjt;
        list<SelRes> ansi,ansj;
        align_segs(i1,maxit,j1,maxjt,ansi,ansj);
        while (maxit!=i2 && maxjt!=j2 && maxit->nam == maxjt->nam) {
            ansi.push_back(*maxit++);
            ansj.push_back(*maxjt++);
        }
        align_segs(maxit,i2,maxjt,j2,ansi,ansj);
        if (ansi.size()==ansj.size() && ansi.size()>3) {
            resi.splice(resi.end(),ansi);
            resj.splice(resj.end(),ansj);
        }
    }

    void RMSD_Utils::seq_align(list<SelRes> &l1,list<SelRes> &l2)
    {
        if (l1.empty() || l2.empty()) return;

        Logger blog;

        list<SelRes> ansi,ansj;
        align_segs(l1.begin(),l1.end(),l2.begin(),l2.end(),ansi,ansj);
        l1.clear();
        l2.clear();
        l1.splice(l1.end(),ansi);
        l2.splice(l2.end(),ansj);

        if (l1.size()!=l2.size()) {
            prf::cerr<<"RMSD_Utils> Error in sequence alignment.\n";
            exit(1);
        }

        blog(log_thres)<<"RMSD_Utils> Sequence alignment: \n";

        for (list<SelRes>::iterator it=l1.begin(),jt=l2.begin();
            it!=l1.end();++it,++jt) {
            blog<<it->nam<<"("<<it->indx_str<<") aligned with "
            <<jt->nam<<"("<<jt->indx_str<<")\n";
        }
    }

    bool RMSD_Utils::seq_match(list<SelRes> &l1,list<SelRes> &l2)
    {
        if (l1.empty()||l2.empty()) return true;

        Logger blog(log_thres);
        if (l1.size()!=l2.size()) {
            blog<<"RMSD_Utils> sequence sizes "<<l1.size()<<" and "<<l2.size()
            <<"\n";
            return false;
        }

        bool matched=true;

        list<SelRes>::iterator it,jt;

        for (it=l1.begin(),jt=l2.begin();it!=l1.end();++it,++jt) {
            if (!(matched=it->nam==jt->nam)) {
                prf::cerr<<"RMSD_Utils> Mismatch at "<<it->nam<<"("<<it->indx_str
                <<") vs "<<jt->nam<<"("<<jt->indx_str<<")\n";
                break;
            }
        }

        return matched;
    }

    bool RMSD_Utils::clean_seq(std::list<SelRes> &tmpres)
    {
        Logger blog(log_thres/2);
        int loc=0;
        while (!(tmpres.empty()||Groups::checkGroup(tmpres.front().nam))) {
            blog<<"While constructing protein sequence, dropped unusable group "
            <<tmpres.front().nam<<". \n";
            loc=1;
            tmpres.pop_front();
        }
        while (!(tmpres.empty()||Groups::checkGroup(tmpres.back().nam))) {
            blog<<"While constructing protein sequence, dropped unusable group "
            <<tmpres.back().nam<<". \n";
            tmpres.pop_back();
        }
        bool cleanseq=true;
        for (list<SelRes>::iterator lr=tmpres.begin();lr!=tmpres.end();++lr) {
            if (!Groups::checkGroup(lr->nam)) {
                prf::cerr<<"The selection contains at least one unknown group "
                <<(lr->nam)<<" at index "<<loc
                <<". The unknown group does not appear at a chain end, and "
                <<"hence, no continuous sequence can be created from this "
                <<"selection.\n";
                cleanseq=false;
                break;
            }
            ++loc;
        }
        return cleanseq;
    }
}
