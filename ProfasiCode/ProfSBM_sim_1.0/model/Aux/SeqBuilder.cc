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

#include "SeqBuilder.hh"
#include "fileutils.hh"

using std::vector;
using std::string;

namespace prf
{
    SeqBuild::SeqBuild() {charmode=false;workstr="";}

    SeqBuild::~SeqBuild() {}

    char SeqBuild::next_char(int &ist)
    {
        char c='?';

        while ((c=workstr[ist])&&isspace(c)) ++ist;

        if (c) ++ist;

        return c;
    }

    string SeqBuild::next_word(int &ist)
    {
        string wd="";
        char c;

        while ((c=workstr[ist])&&isspace(c)) ist++;

        while ((c=workstr[ist])&&!(isspace(c)||c=='*')) {
            wd.push_back(c);ist++;
        }

        return wd;
    }

    OneLetterCode SeqBuild::next_code(int &ist)
    {
        OneLetterCode cd=NONE;

        if (charmode) {
            char c=next_char(ist);

            if (!c) return cd;

            cd=Groups::mapChar2OLC(c);

            if (c=='*') {
                charmode=false;
                return next_code(ist);
            }

            if (cd==NONE) prf::cerr<<"Unknown group "<<c<<"\n";

            return cd;
        } else {
            string s=next_word(ist);

            if (!s.empty()) {
                Groups::checkGroup(s);
                cd=Groups::map2OLC(s);
            }

            char c=next_char(ist);

            if (c=='*') charmode=true;

            if (c && c!='*') --ist;

            return cd;
        }
    }

    void SeqBuild::parse(string inpstr,vector<OneLetterCode> &sq)
    {
        sq.clear();
        int i=0,j=0;
        char c;
        inpstr=prf_utils::trim_str(inpstr);

        if (inpstr[i]=='<') ++i;
        j=inpstr.size()-1;
        if (j<i) return;
        if (inpstr[j]=='>') --j;
        if (j<i) return;
        while (isspace(c=inpstr[j]) || c=='*') --j;

        if (j<i) return;

        workstr=inpstr.substr(i,j-=(i-1));

        for (i=0;i<j;) {
            OneLetterCode cd=next_code(i);

            if (cd!=NONE) sq.push_back(cd);

            if (!workstr[i]) break;
        }

        charmode=false;
    }

    std::string SeqBuild::make_string(const std::vector<OneLetterCode> &sq,
                                char begn, char nd)
    {
        std::string ans;
        bool prevmode=charmode;
        word_mode();
        for (size_t i=0;i<sq.size();++i) {
            if (Groups::grp[sq[i]].isAA() && Groups::grp[sq[i]].isNat()) {
                if (!charmode) {
                    letter_mode();
                    ans+=" * ";
                }
                ans+=Groups::mapOLC2Char(sq[i]) ;
            } else {
                if (charmode) {
                    word_mode();
                    ans+=" * ";
                }
                ans+=Groups::mapOLC2TLC(sq[i]);
            }
        }
        if (charmode!=prevmode) ans+=" * ";
        charmode=prevmode;
        return begn+ans+nd;
    }
}
