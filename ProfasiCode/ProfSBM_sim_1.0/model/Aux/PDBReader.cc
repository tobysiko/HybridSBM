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

#include "PDBReader.hh"
#include <fstream>
#include <algorithm>
#include <deque>
#include <set>
using std::ifstream;
using std::string;
using std::list;
using std::deque;
using std::vector;

using namespace prf_utils;

namespace prf
{

    using namespace prf_pdb_vars;
    PDB_model::PDB_model()
    {
        num_chains=0;
        chain_name.clear();
        res_map.clear();
    }

    PDB_model::~PDB_model() {}

    PDB_model::PDB_model(const PDB_model &pm)
    {
        num_chains=0;
        chain_name=pm.chain_name;
        res_map=pm.res_map;
    }

    PDB_model & PDB_model::operator=(const PDB_model &pm)
    {
        if (this!=&pm) {
            num_chains=0;
            chain_name=pm.chain_name;
            res_map=pm.res_map;
        }

        return *this;
    }

    int PDB_model::index_of_res(string chid, string resnum)
    {
        int ich=index_of(chid);

        for (size_t j=0;j<res_map[ich].size();++j) {
            if (res_map[ich][j]==resnum) return j;
        }

        return -1;
    }

    int PDB_model::index_of(string chid)
    {
        for (size_t i=0;i<chain_name.size();++i) {
            if (chid==chain_name[i]) return i;
        }

        return 0;
    }

    void PDB_model::print_summary(Output &op)
    {
        op<<"Number of chains = "<<num_chains<<"\n";
        op<<"Chain name, residue index, file residue index, matrix lines\n";

        for (int i=0;i<num_chains;++i) {
            op<<chain_name[i]<<"...\n";

            for (size_t j=0;j<res_map[i].size();++j) {
                op<<j<<"  "<<res_map[i][j]<<"  "
                <<res_ln0[i][j]<<"   "<<res_lnn[i][j]<<"\n";
            }
        }
    }

    PDBReader::PDBReader() : PopBase()
    {
        matrix_read=false;
        curmodel=0;
        nto=0;
    }

    PDBReader::PDBReader(string fl) : PopBase()
    {
        matrix_read=false;
        curmodel=0;
        nto=0;
        set_file(fl);
    }

    PDBReader::~PDBReader() {}

    void PDBReader::set_file(string fl)
    {
        myfile=fl;
    }

    void PDBReader::set_model(int i)
    {
        PopBase::set_model(i);
        nc=models[curmodel].num_chains;
    }

    string PDBReader::str_index(int ires, int ic)
    {
        return models[curmodel].res_map[ic][ires];
    }

    string PDBReader::grp_name(int i, int ich)
    {
        string tlc =read_word(models[curmodel].res_ln0[ich][i],
                              res_label_column,res_label_width);
        return tlc;
    }

    int PDBReader::sequence_of_chain(int ich,list<string> &seq)
    {
        seq.clear();
        int res0=0,resn=models[curmodel].res_map.size();

        for (int i=res0;i<resn;++i) {
            seq.push_back(grp_name(i,ich));
        }

        return (int) seq.size();
    }

    string PDBReader::read_word(int rwo, int stcol, int len)
    {
        string wd;

        for (int j=0;j<len;++j) wd+=data[rwo][stcol+j];

        return wd;
    }

    void PDBReader::write_word(int rwo, int stcol, std::string txt)
    {
        for (size_t j=0;j<txt.size();++j) data[rwo][stcol+j]=txt[j];
    }

    int PDBReader::clear()
    {
        return 1;
    }

    int PDBReader::read_matrix()
    {
        if (matrix_read) return 2;

        Logger blog;

        nto=0;

        if (TestFile_r(myfile)!=0) {
            ifstream fin(myfile.c_str());
            string kwd,ln;

            while (getline(fin,ln)) ++nto;

            fin.close();

            blog(50)<<"PDBReader> counted "<<nto<<" lines.\n";

            data.allocate(nto,80);

            for (int i=0;i<nto;++i) {
                for (int j=0;j<80;++j) {
                    data[i][j]=' ';
                }
            }

            blog(50)<<"PDBReader> allocated and filled matrix with space\n";

            ifstream fin1(myfile.c_str());
            int curln=0,curcl=0;

            while (getline(fin1,ln)) {
                curcl=0;

                for (string::iterator it=ln.begin();it!=ln.end()&&curcl<80;++it) {
                    data[curln][curcl++]=*it;
                }

                ++curln;
            }

            fin1.close();

            gather_model_info();
            blog<<"PDBReader> Model information gathering finished.\n";
            check_file();
            matrix_read=true;
            return 1;
        }

        matrix_read=false;

        return 0;
    }

    void PDBReader::delete_matrix()
    {
        if (matrix_read) {
            data.allocate(2,2);
            models.clear();
            matrix_read=false;
        }
    }

    void PDBReader::gather_model_info()
    {
        // count models and mark model boundaries
        Logger blog;
        nmodels=0;
        vector<int> mod_begin,mod_end,mod_cl,ch_begin,ch_end;
        string kwd;
        mod_begin.clear();mod_end.clear();
        nmodels=0;

        for (int i=0;i<nto;++i) {
            kwd=read_word(i,0,6);

            if (kwd=="MODEL ") {
                mod_begin.push_back(i);
                ++nmodels;
                mod_end.push_back(-1);
                mod_cl.push_back(0);
            }

            if (kwd=="ENDMDL") {
                mod_end[nmodels-1]=i;
                mod_cl[nmodels-1]=1;
            }
        }

        blog(30)<<"PDBReader> Number of \"MODEL\" statements in "<<myfile

        <<" = "<<nmodels<<"\n";

        for (int i=0;i<nmodels;++i) {
            if (mod_cl[i]!=1) {
                blog(30)<<"PDBReader> File "<<myfile<<": model "<<i
                <<", closure error.\n";

                if (i!=(nmodels-1)) mod_end[i]=mod_begin[i+1]-1;
                else mod_end[i]=nto;
            }
        }

        if (nmodels==0) {
            blog(20)<<"PDBReader> File "<<myfile<<" does not have MODEL tags. "
                    <<"This means, the file contains only one model.\n";
            nmodels=1;
            mod_begin.push_back(0);
            mod_end.push_back(nto);
            mod_cl.push_back(1);
        } else {
            blog(5)<<"PDBReader> File "<<myfile<<" contains "<<nmodels
                    <<" models\n";
        }

        models.resize(nmodels);

        blog(50)<<"PDBReader> File "<<myfile<<":\n";

        for (int i=0;i<nmodels;++i) {
            blog<<"Model "<<i<<" begins at "<<mod_begin[i]
                    <<" ends at "<<mod_end[i]<<"\n";
        }

        for (int i=0;i<nmodels;++i) {
            // for each model, find how many chains there are and their labels
            int nch=0;
            string chnlbl="",chnlblp="",resnum="",resnump="";
            ch_begin.clear();ch_end.clear();

            for (int j=mod_begin[i];j<mod_end[i];++j) {
                kwd=read_word(j,0,6);

                if (kwd!="ATOM  " && kwd!="HETATM") continue;

                chnlbl=read_word(j,chain_label_column,chain_label_width);

                if (chnlbl==" ") chnlbl="A";

                if (chnlbl!=chnlblp) {
                    ++nch;
                    models[i].chain_name.push_back(chnlbl);
                    ch_begin.push_back(j);
                    ch_end.push_back(mod_end[i]);

                    if (nch>1) ch_end[nch-2]=j;

                    chnlblp=chnlbl;

                    blog(50)<<"PDBReader> Read new chain label "<<chnlbl<<"\n";
                }
            }
            ch_end[nch-1]=mod_end[i];
            for (int k=0;k<nch;++k) {
                int prevline=0;
                chnlbl=models[i].chain_name[k];
                do {
                    prevline=ch_end[k]-1;
                    if (prevline<ch_begin[k]) break;
                    kwd=read_word(prevline,0,6);
                    chnlblp=read_word(ch_begin[k],chain_label_column,chain_label_width);
                    if (chnlblp==" ") chnlblp="A";
                    if ((kwd!="ATOM  " && kwd!="HETATM") or chnlblp!=chnlbl) ch_end[k]=prevline;
                    else break;
                } while (1);
            }
            models[i].num_chains=nch;

            models[i].res_map.resize(nch);
            models[i].res_ln0.resize(nch);
            models[i].res_lnn.resize(nch);

            for (int ich=0;ich<nch;++ich) {
                resnump="";
                for (int k=ch_begin[ich];k<ch_end[ich];++k) {
                    kwd=read_word(k,0,6);

                    if (kwd!="ATOM  " && kwd!="HETATM") continue;

                    resnum=read_word(k,res_num_column,res_num_width);

                    if (resnum!=resnump) {
                        models[i].res_map[ich].push_back(resnum);
                        resnump=resnum;
                        models[i].res_ln0[ich].push_back(k);
                    }
                }

                models[i].res_lnn[ich].resize(models[i].res_ln0[ich].size(),0);

                for (size_t k=0;k<models[i].res_lnn[ich].size()-1;++k) {
                    models[i].res_lnn[ich][k]=models[i].res_ln0[ich][k+1];
                }

                models[i].res_lnn[ich][models[i].res_lnn[ich].size()-1]
                =ch_end[ich];

            }
        }

        set_model(0);
    }

    string PDBReader::chain_name(int i) const
    {
        return models[curmodel].chain_name[i];
    }

    int PDBReader::num_grp(int ich) const
    {
        return (int) models[curmodel].res_map[ich].size();
    }

    int PDBReader::index_of_grp(string ires, int ich)
    {
        size_t i=0;

        while (ires.size()<3) ires=" "+ires;

        for (;i<models[curmodel].res_map[ich].size();++i)
            if (models[curmodel].res_map[ich][i]==ires) break;

        return (int) i;
    }

    Vector3 PDBReader::mk_vec(int lno)
    {
        return Vector3(atof(read_word(lno,crd_x_col,crd_col_width).c_str()),
                       atof(read_word(lno,crd_y_col,crd_col_width).c_str()),
                       atof(read_word(lno,crd_z_col,crd_col_width).c_str()));
    }

    void PDBReader::mk_rec(int lno, AtomRecord &ds)
    {
        char line[81];

        line[80]=0;

        for (int k=0;k<80;++k) line[k]=data[lno][k];

        ds.set_fields(string(line), &aldict);

        ds.descriptor().int_label=lno;
    }

    int PDBReader::add_descriptors(SelRes &sel,
                                   std::list<AtomDescriptor> &des)
    {
        int ich=sel.chn;
// set_model(sel.mdl);
        int fstline=models[curmodel].res_ln0[ich][sel.indx_nat];
        int lstline=models[curmodel].res_lnn[ich][sel.indx_nat];

        for (int j=fstline;j<lstline;++j) {
            AtomRecord ard;
            mk_rec(j,ard);
            ard.descriptor().iresrel=sel.indx_nat;
            des.push_back(ard.descriptor());
        }

        return lstline-fstline;
    }

    int PDBReader::add_records(SelRes &sel, std::list<AtomRecord> &des)
    {
        int ich=sel.chn;
// set_model(sel.mdl);
        int fstline=models[curmodel].res_ln0[ich][sel.indx_nat];
        int lstline=models[curmodel].res_lnn[ich][sel.indx_nat];

        for (int j=fstline;j<lstline;++j) {
            AtomRecord ard;
            mk_rec(j,ard);
            ard.descriptor().iresrel=sel.indx_nat;
            des.push_back(ard);
        }

        return lstline-fstline;
    }

    int PDBReader::descriptors(list<SelRes> &slc,
                               list<AtomDescriptor> &des)
    {
        int nadded=0;
        des.clear();

        for (list<SelRes>::iterator i=slc.begin();i!=slc.end();++i) {
            nadded+=add_descriptors(*i,des);
        }

        return nadded;
    }

    int PDBReader::records(list<SelRes> &slc, list<AtomRecord> &des)
    {
        int nadded=0;
        des.clear();

        for (list<SelRes>::iterator i=slc.begin();i!=slc.end();++i) {
            nadded+=add_records(*i,des);
        }

        return nadded;
    }

    int PDBReader::export_shape(vector<int> &vct, Shape &shp)
    {
        shp.Resize(vct.size());

        for (size_t i=0;i<vct.size();++i) {
            shp.Point(i,mk_vec(vct[i]));
        }

        return (int) vct.size();
    }

    int PDBReader::export_records(vector<int> &vct,
                                  list<AtomRecord> &shp)
    {
        shp.clear();

        for (size_t i=0;i<vct.size();++i) {
            AtomRecord ard;
            mk_rec(vct[i],ard);
            shp.push_back(ard);
        }

        return 1;
    }

    int PDBReader::export_descriptors(vector<int> &vct,
                                      list<AtomDescriptor> &shp)
    {
        shp.clear();

        for (size_t i=0;i<vct.size();++i) {
            AtomRecord ard;
            mk_rec(vct[i],ard);
            shp.push_back(ard.descriptor());
        }

        return 1;
    }

    bool PDBReader::check_file()
    {
        bool allok=true;
        allok = allok && check_unknown_labels();
        return allok;
    }

    bool PDBReader::check_unknown_labels()
    {
        bool labels_ok=false;
        Logger blog(20);
        std::vector<std::set<std::string> > bad_labels(Groups::max_olc);
        std::vector<std::set<std::string> > hlabels(Groups::max_olc);
        std::map<std::string,int> ignored;
        for (int ich=0;ich<nc;++ich) {
            for (size_t j=0; j<models[curmodel].res_map[ich].size();++j) {
                int fstline=models[curmodel].res_ln0[ich][j];
                int lstline=models[curmodel].res_lnn[ich][j];
                string reslabel=grp_name(j,ich);
                OneLetterCode curgrp=Groups::map2OLC(reslabel);
                if (curgrp==NONE) {
                    ignored[reslabel]=ignored[reslabel]+1;
                    continue;
                }
                for (int k=fstline;k<lstline;++k) {
                    std::string kwd=read_word(k,0, 6);

                    if (kwd!="ATOM  " and kwd!="HETATM") continue;

                    string alabel=read_word(k,atom_label_column, atom_label_width);

                    if (!aldict.interpret(curgrp,alabel)) {
                        bad_labels[curgrp].insert(alabel);
                    }
                }

                if (!bad_labels[curgrp].empty()) {
                    blog<<"Found "<<bad_labels[curgrp].size()<<" unfamiliar labels for "
                    <<reslabel<<"\n{";
                    for (std::set<std::string>::iterator k=bad_labels[curgrp].begin();
                    k!=bad_labels[curgrp].end();++k) {
                        blog<<"\""<<(*k)<<"\", ";
                    }
                    blog<<"}\nI will now try to figure out what they might be...\n";

                    for (int k=fstline;k<=lstline;++k) {
                        std::string kwd=read_word(k,0, 6);

                        if (kwd!="ATOM  " and kwd!="HETATM") continue;

                        string alabel=read_word(k,atom_label_column, atom_label_width);

                        if (aldict.is_hydrogen_label(alabel)) {
                            hlabels[curgrp].insert(alabel);
                        }
                    }

                    if (!hlabels[curgrp].empty()) aldict.apply_heuristics(curgrp, hlabels[curgrp]);

                    std::set<std::string>::iterator it,jt=bad_labels[curgrp].begin();

                    while (jt!=bad_labels[curgrp].end()) {
                        it=jt++;
                        std::string tmpstr=*it;
                        if (aldict.interpret(curgrp,tmpstr)) bad_labels[curgrp].erase(it);
                    }

                    blog<<"Remaining number of unfamiliar labels for "
                            <<reslabel<<" = "<<bad_labels[curgrp].size()<<"\n";
                }
                if (curgrp==R) {
                    std::map<std::string, Vector3> crd;
                    std::map<std::string, int> lineof;
                    for (int k=fstline;k<lstline;++k) {
                        std::string kwd=read_word(k,0, 6);

                        if (kwd!="ATOM  " and kwd!="HETATM") continue;

                        string alabel=read_word(k,atom_label_column,
                                                atom_label_width);

                        if (aldict.interpret(curgrp,alabel)) {
                            crd[alabel]=mk_vec(k);
                            lineof[alabel]=k;
                        }
                    }
                    if (crd.find(" CD ")!=crd.end() and
                        crd.find(" NE ")!=crd.end() and
                        crd.find(" CZ ")!=crd.end() and
                        crd.find(" NH2")!=crd.end()) {

                        Vector3 v1,v2,v3;
                        v1=crd[" NE "]-crd[" CD "];
                        v2=crd[" CZ "]-crd[" NE "];
                        v3=crd[" NH2"]-crd[" CZ "];
                        if (v1.dot(v3)>0) {
                            write_word(lineof[" NH2"],atom_label_column," NH1");
                            if (lineof.find(" NH1")!=lineof.end()) {
                                write_word(lineof[" NH1"],atom_label_column,
                                           " NH2");
                            }
                            if (lineof.find("1HH1")!=lineof.end()) {
                                write_word(lineof["1HH1"],atom_label_column,
                                           "1HH2");
                            }
                            if (lineof.find("2HH1")!=lineof.end()) {
                                write_word(lineof["2HH1"],atom_label_column,
                                           "2HH2");
                            }
                            if (lineof.find("1HH2")!=lineof.end()) {
                                write_word(lineof["1HH2"],atom_label_column,
                                           "1HH1");
                            }
                            if (lineof.find("2HH2")!=lineof.end()) {
                                write_word(lineof["2HH2"],atom_label_column,
                                           "2HH1");
                            }
                        }
                    }
                    //Ok, this is hideous. This class needs a redesign.
                }
            }
        }

        if (!ignored.empty()) {
            prf::cerr<<"\nWhile inspecting atom labels, the following groups were ignored.\n...\n";
            for (std::map<std::string,int>::iterator it=ignored.begin();
                    it!=ignored.end();++it) {
                prf::cerr<<it->second<<" instances of group "<<it->first<<"\n";
            }
            prf::cerr<<"...\nThese groups have not been modelled in PROFASI.\n\n";
        }

        int nunknown=0;
        for (int i=0;i<Groups::max_olc;++i) nunknown+=bad_labels[i].size();

        if (nunknown!=0) {
            prf::cerr<<"PDBReader> Found "<<nunknown<<" unknown atom labels. \n";
            prf::cerr<<"Use the pdb_energies program in the interactive mode to "
                <<"train the label interpreter.\n";
        } else labels_ok=true;

        return labels_ok;
    }
}
