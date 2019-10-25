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

#include "AtomLabelDictionary.hh"
#include <cstdlib>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <errno.h>
#include <limits.h>
#include "profasi_io.hh"
#include "fileutils.hh"
#include "../Elements/GroupLib.hh"

using namespace prf_utils;

using namespace prf::Groups;

namespace prf
{
    int flush_cin()
    {
        int nxt=std::cin.peek();

        while (nxt=='\n') {
            std::cin.get();
            nxt=std::cin.peek();
        }

        return nxt;
    }

    int get_option()
    {
        flush_cin();
        std::string tmpstr;
        getline(std::cin,tmpstr);
        errno=0;
        long vl=strtol(tmpstr.c_str(),NULL,10);

        if (errno!=0) vl=INT_MIN;

        return vl;
    }

    AtomLabelDictionary::AtomLabelDictionary()
    {
        char *homeDir = getenv("HOME");
        std::string filename = (homeDir == 0) ?
                              std::string("./.ProFASi_dictionary"):
                              homeDir + std::string("/.ProFASi_dictionary");

        if (prf_utils::STestFile_r(filename.c_str())) {
            std::ifstream fin(filename.c_str());
            std::string line;

            while (getline(fin,line)) {
                if (line.empty()) continue;

                std::vector<std::string> parts;
                split_str(line,'/',parts);

                if (parts.size()<3) continue;

                prf::OneLetterCode olc=map2OLC(parts[0]);

                if (olc==NONE) continue;

                if (Groups::grp[olc].valid_label(parts[2])) {
                    lblmap[olc][parts[1]]=parts[2];
                }
            }

            fin.close();
        }

        termlbl["1H"]="1H  ";
        termlbl["2H"]="2H  ";
        termlbl["3H"]="3H  ";
        termlbl["H1"]="1H  ";
        termlbl["H2"]="2H  ";
        termlbl["H3"]="3H  ";
        termlbl["OXT"]=" OXT";
    }

    AtomLabelDictionary::~AtomLabelDictionary()
    {
        if (!lblcache.empty()) {
            prf::cout<<"\n\nThere are some new instructions on atom labels ...\n";

            for (std::map<OneLetterCode, std::map<std::string,std::string> >::iterator
                 it=lblcache.begin(); it!=lblcache.end(); ++it) {
                if (!(it->second.empty())) {
                    for (std::map<std::string, std::string>::iterator
                         jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                        prf::cout<<"For "<<mapOLC2Name(it->first)
                                 <<", when I see \""<<jt->first
                                 <<"\", I should read \""<<jt->second<<"\". \n";
                    }
                }
            }

            char save='y';
            prf::cout<<"\n\nDo you want to remember this for the future use ? (y/n):";
            std::cin>>save;

            if (save=='y') {
                std::string filename=getenv("HOME")+std::string("/.ProFASi_dictionary");
                Output dct(filename.c_str(),"a");

                for (std::map<OneLetterCode, std::map<std::string,std::string> >::iterator
                     it=lblcache.begin(); it!=lblcache.end(); ++it) {
                    if (!(it->second.empty())) {
                        for (std::map<std::string, std::string>::iterator
                             jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                            dct<<mapOLC2Name(it->first)<<"/"<<jt->first<<"/"
                            <<jt->second<<"/\n";
                        }
                    }
                }

                dct.close();
                prf::cout<<"Mappings saved in file "<<filename<<" for future reference.\n";
                prf::cout.flush();
            }
        }
    }

    bool AtomLabelDictionary::interpret(std::string tlc, std::string &alabel)
    {
        return interpret(map2OLC(tlc), alabel);
    }

    bool AtomLabelDictionary::interpret(OneLetterCode olc, std::string &alabel)
    {
        std::string dummy;

        if (!lblmap[olc][alabel].empty()) {
            alabel=lblmap[olc][alabel];
            return true;
        }

        //User specified mappings take precedence over program defaults
        if (Groups::grp[olc].valid_label(alabel)) return true;

        dummy=rm_space(alabel);

        if (!termlbl[dummy].empty()) {
            alabel=termlbl[dummy];
            return true;
        }

        return false;
    }

    int AtomLabelDictionary::learn(std::vector<std::set<std::string> > &new_labels)
    {
        std::vector<prf::OneLetterCode> relevant_res;
        std::vector<int> num_unknown;
        int nnew=0; //The number of interpretations changed interactively

        for (size_t i=0; i<new_labels.size(); ++i) {
            if (!new_labels[i].empty()) {
                relevant_res.push_back((prf::OneLetterCode) i);
                num_unknown.push_back(0);
            }
        }

        if (relevant_res.empty()) return true;

        if (relevant_res.size()==1)
            return (learn_res(relevant_res[0], new_labels[relevant_res[0]])==0);

        prf::cout<<"The learn function was called with labels for the following groups...\n";
        bool finished=false;

        do {
            for (size_t ir=0; ir<relevant_res.size(); ++ir) {
                int nunk=0;

                for (std::set<std::string>::iterator li=new_labels[relevant_res[ir]].begin();
                     li!=new_labels[relevant_res[ir]].end(); ++li) {
                    std::string tmpstr=*li;

                    if (!interpret(relevant_res[ir],tmpstr)) ++nunk;
                }

                num_unknown[ir]=nunk;
                prf::cout<<ir<<": "<<Groups::mapOLC2Name(relevant_res[ir])
                         <<"("<<Groups::mapOLC2TLC(relevant_res[ir])<<") ";

                if (num_unknown[ir]!=0) prf::cout<<"["<<num_unknown[ir]
                                                     <<" unknown atom labels.]\n";
                else prf::cout<<"\n";
            }

            prf::cout<<"Choose one to view/edit labels (-1, to quit): ";
            int jres=get_option();

            if (jres<-1 or jres>=(int)relevant_res.size()) continue;

            if (jres==-1) finished=true;
            else nnew+=learn_res(relevant_res[jres], new_labels[relevant_res[jres]]);
        } while (not finished);

        return nnew;
    }

    int AtomLabelDictionary::learn_res(prf::OneLetterCode olc,
                                       std::set<std::string> &lbls)
    {
        std::string resn=mapOLC2Name(olc);
        int nunknown=0,ic=0,nnew=0;
        std::vector<std::string> known=Groups::grp[olc].valid_labels();
        std::vector<std::string> lblsv(lbls.size());

        for (std::set<std::string>::iterator i=lbls.begin(); i!=lbls.end(); ++i) {
            lblsv[ic++]=*i;
        }

        do {
            nunknown=list_interpretations(olc,lblsv);
            prf::cout<<"\n";
            prf::cout<<"If you wish to change the way a label is interpreted, enter it's "
                     <<"number. (Enter -1 to quit this menu):";
            int choice1=get_option();

            if (choice1==-1) break;
            else if (choice1<0 or choice1>= (int) lblsv.size()) continue;

            prf::cout<<"The following are the list of atom labels for "<<resn
                     <<", recognized by PROFASI ...\n";

            for (size_t i=0; i<known.size(); ++i) {
                prf::cout<<i<<":  \""<<known[i]<<"\"\t\t";

                if ((i+1)%5==0) prf::cout<<"\n";
            }

            prf::cout<<"\nTo map \""<<lblsv[choice1]
                     <<"\" to one of the above, enter its number. \n"
                     <<"(Enter -1 to quit, -2 for further mapping options):";
            int choice2=get_option();

            if (choice2==-2) {
                prf::cout<<"You can choose to map "<<lblsv[choice1]<<" to any 4 letter "
                         <<"string of your choice. Spaces in your string should be entered as "
                         <<"underscores ('_'):\n";
                std::string intp, newlabel;
                prf::cout<<"Enter your interpretation for \""<<lblsv[choice1];
                prf::cout<<"\" (empty string leaves it unchanged):";
                flush_cin();
                getline(std::cin,intp);
                intp=rm_space(intp);

                if (!intp.empty()) {
                    newlabel=std::string(intp,0,4);

                    for (int jj=0; jj<4; ++jj) if (newlabel[jj]=='_') newlabel[jj]=' ';

                    prf::cout<<"Is this alias specific for residue "<<resn
                             <<"(0) or generic (1) ? ";
                    int choice3=get_option();

                    if (choice3==0) {
                        nnew+=make_alias(olc,lblsv[choice1],intp);
                    } else {
                        size_t oldsz=termlbl.size();
                        termlbl[rm_space(lblsv[choice1])]=newlabel;
                        nnew+=termlbl.size()==oldsz?0:1;
                    }
                }
            } else if (choice2>=0 and choice2<(int)known.size()) {
                nnew+=make_alias(olc,lblsv[choice1],known[choice2]);
                lblcache[olc][lblsv[choice1]]=known[choice2];
            }
        } while (1);

        return nnew;
    }
    int AtomLabelDictionary::list_interpretations(prf::OneLetterCode olc,
            std::vector<std::string> &lbls)
    {
        int nunknown=0;
        prf::cout<<"Residue:"<<Groups::mapOLC2Name(olc)
                 <<", labels with their interpretations (in paranthesis)\n";

        for (size_t i=0; i<lbls.size(); ++i) {
            std::string interp=lbls[i];
            prf::cout<<i<<". \""<<lbls[i]<<"\"";

            if (interpret(olc,interp)) {
                prf::cout<<"(\""<<interp<<"\")\t";
            } else {
                prf::cout<<"(\"????\")\t";
                ++nunknown;
            }

            if ((i+1)%5==0) prf::cout<<"\n";
        }

        return nunknown;
    }
    bool AtomLabelDictionary::is_hydrogen_label(std::string alabel)
    {
        return (alabel[0]=='H' or
                ((isdigit(alabel[0]) or alabel[0]==' ') and alabel[1]=='H'));
    }

    void AtomLabelDictionary::apply_heuristics(prf::OneLetterCode olc,
            std::set<std::string> &hlbls)
    {
        std::vector<std::string> mylbls=Groups::grp[olc].valid_labels();
        std::list<std::string> myhlbls,myhclass, hclass;

        for (size_t i=0; i<mylbls.size(); ++i) {
            if (is_hydrogen_label(mylbls[i])) {
                myhlbls.push_back(mylbls[i]);
            }
        }

        while (!myhlbls.empty()) {
            myhclass.clear();
            hclass.clear();

            do {
                myhclass.splice(myhclass.end(),myhlbls,myhlbls.begin());
            } while ((!myhlbls.empty()) and
                     (class_of(myhlbls.front())==class_of(myhclass.back())));

            std::string curclass=class_of(myhclass.front());
            std::set<std::string>::iterator i1,i2;
            i1=i2=hlbls.begin();

            while (i2!=hlbls.end()) {
                i1=i2++;

                if (!curclass.empty()) {
                    if (i1->find(curclass)<i1->size()) {
                        hclass.push_back(*i1);
                        hlbls.erase(i1);
                    }
                } else {
                    int n_al=0, n_num=0, n_h=0;

                    for (size_t jl=0; jl<i1->size(); ++jl) {
                        if ((*i1)[jl]=='H') ++n_h;

                        if (isalpha((*i1)[jl])) ++n_al;

                        if (isdigit((*i1)[jl])) ++n_num;
                    }

                    if (n_al==1 && n_h==1) {
                        hclass.push_back(*i1);
                        hlbls.erase(i1);
                    }
                }
            }

            std::list<std::string>::iterator ii1,ii2;
            ii1=hclass.begin();
            ii2=myhclass.begin();
            int smallerlist=(int)hclass.size();

            if (myhclass.size()<hclass.size()) smallerlist=myhclass.size();

            while (smallerlist-->0) make_alias(olc,*ii1++,*ii2++);
        }
    }

    std::string AtomLabelDictionary::class_of(std::string anhlbl)
    {
        std::string tmp = anhlbl.substr(2);

        while (isspace(tmp[tmp.size()-1])) tmp.resize(tmp.size()-1);

        return tmp;
    }

    int AtomLabelDictionary::make_alias(prf::OneLetterCode olc,
                                        std::string s1, std::string s2)
    {
        if (s1==s2) return 0;

        std::string tmp=s1;
        interpret(olc,tmp);
        Logger blog(20);
        blog<<"In residue "<<Groups::mapOLC2Name(olc)<<", ";

        if (tmp!=s2) {
            blog<<"creating alias \""<<s1<<"\" for \""<<s2<<"\"\n";
            lblmap[olc][s1]=s2;
            return 1;
        } else {
            blog<<"alias \""<<s1<<"\" for \""<<s2<<"\" was not created, as "<<s1
            <<" is already interpreted like that.\n";
            return 0;
        }
    }
}
