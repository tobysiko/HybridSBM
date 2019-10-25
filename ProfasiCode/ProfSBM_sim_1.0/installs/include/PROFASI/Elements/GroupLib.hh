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

#ifndef GroupLib_HH
#define GroupLib_HH
#include "GroupProps.hh"
#include "../Aux/profasi_io.hh"

//Automatically generated group library header file

namespace prf
{

    typedef enum { NONE, G, A, V, L, I, S, T, C, M, P, D, N, E, Q, K, R, H, F, Y, W, DPR, ACE, NH2, SUC, NME, VOIDEG } OneLetterCode;

    namespace Groups
    {
        const int max_olc=27;
        extern GroupProps grp[max_olc];
        extern char charcode[max_olc];
        extern std::map<std::string,OneLetterCode> olcof;
        bool checkGroup(std::string gg);
        void initGroups() ;


        inline OneLetterCode mapTLC2OLC(std::string tlc) {return olcof[tlc];}

        inline OneLetterCode map2OLC(std::string nm) {return olcof[nm];}

        inline bool isAA(std::string st) {return grp[map2OLC(st)].isAA();}

        inline bool isEG(std::string st){return grp[map2OLC(st)].isEG();}

        inline std::string egname(std::string st)
        {
            if (grp[map2OLC(st)].isEG()) return grp[map2OLC(st)].CommonName();
            else return std::string("null");
        }

        inline Output & operator<<(Output & os, OneLetterCode olc) {return os<<((int) olc);}

        inline char mapOLC2Char(OneLetterCode olc) { return charcode[olc];}

        inline std::string mapOLC2Name(OneLetterCode olc) { return grp[olc].CommonName();}

        inline std::string mapOLC2TLC(OneLetterCode olc) { return grp[olc].TLC();}

        inline std::string SCALabel(OneLetterCode olc,int i) { return grp[olc].label(grp[olc].index(" CB ")+i);}

        inline OneLetterCode mapChar2OLC(char c) {return olcof[std::string(1,c)];}

    }
}

#endif

