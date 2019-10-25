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

#include "GroupLib.hh"

namespace prf
{

    namespace Groups {

        GroupProps grp[max_olc];
        char charcode[max_olc];
        std::map<std::string,OneLetterCode> olcof;
        bool checkGroup(std::string gg)
        {
            if (map2OLC(gg)==NONE) {
                prf::cerr<<"Unknown group "<<gg<<"\n";
                return false;
            }

            return true;
        }

        void initGroups()
        {
            static bool already_initialized=false;

            if (already_initialized) return;

            //Ghost
            grp[NONE].CommonName("Ghost");

            grp[NONE].init(std::string("XXXX"));

            grp[NONE].set_type("unnatural things");

            grp[NONE].TLC("NON");

            charcode[NONE]='X';

            olcof["NON"]=NONE;

            olcof["Ghost"]=NONE;

            //Glycine
            grp[G].CommonName("Glycine");

            grp[G].init(std::string(" N  , H  , CA ,1HA ,2HA , C  , O  "));

            grp[G].set_type("natural aminoacid");

            grp[G].TLC("GLY");

            charcode[G]='G';

            olcof["GLY"]=G;

            olcof["Glycine"]=G;

            olcof["G"]=G;

            //Alanine
            grp[A].CommonName("Alanine");

            grp[A].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB ,3HB ,") +
                        std::string(" C  , O  "));

            grp[A].set_type("natural hydrophobic aminoacid");

            grp[A].TLC("ALA");

            charcode[A]='A';

            olcof["ALA"]=A;

            olcof["Alanine"]=A;

            olcof["A"]=A;

            grp[A].torsion_dof(0,2,4,5);

            //Valine
            grp[V].CommonName("Valine");

            grp[V].init(std::string(" N  , H  , CA , HA , CB , HB , CG1, CG2,") +
                        std::string("1HG1,2HG1,3HG1,1HG2,2HG2,3HG2, C  , O  "));

            grp[V].set_type("natural hydrophobic aminoacid");

            grp[V].TLC("VAL");

            charcode[V]='V';

            olcof["VAL"]=V;

            olcof["Valine"]=V;

            olcof["V"]=V;

            grp[V].torsion_dof(0,2,4,6);

            grp[V].torsion_dof(2,4,6,8);

            grp[V].torsion_dof(2,4,7,11);

            //Leucine
            grp[L].CommonName("Leucine");

            grp[L].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" HG , CD1, CD2,1HD1,2HD1,3HD1,1HD2,2HD2,") +
                        std::string("3HD2, C  , O  "));

            grp[L].set_type("natural hydrophobic aminoacid");

            grp[L].TLC("LEU");

            charcode[L]='L';

            olcof["LEU"]=L;

            olcof["Leucine"]=L;

            olcof["L"]=L;

            grp[L].torsion_dof(0,2,4,7);

            grp[L].torsion_dof(2,4,7,9);

            grp[L].torsion_dof(4,7,9,11);

            grp[L].torsion_dof(4,7,10,14);

            //Isoleucine
            grp[I].CommonName("Isoleucine");

            grp[I].init(std::string(" N  , H  , CA , HA , CB , HB , CG1, CG2,") +
                        std::string("1HG2,2HG2,3HG2,1HG1,2HG1, CD1,1HD1,2HD1,") +
                        std::string("3HD1, C  , O  "));

            grp[I].set_type("natural hydrophobic aminoacid");

            grp[I].TLC("ILE");

            charcode[I]='I';

            olcof["ILE"]=I;

            olcof["Isoleucine"]=I;

            olcof["I"]=I;

            grp[I].torsion_dof(0,2,4,6);

            grp[I].torsion_dof(2,4,6,13);

            grp[I].torsion_dof(2,4,7,8);

            grp[I].torsion_dof(4,6,13,14);

            //Serine
            grp[S].CommonName("Serine");

            grp[S].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , OG ,") +
                        std::string(" HG , C  , O  "));

            grp[S].set_type("natural aminoacid");

            grp[S].TLC("SER");

            charcode[S]='S';

            olcof["SER"]=S;

            olcof["Serine"]=S;

            olcof["S"]=S;

            grp[S].torsion_dof(0,2,4,7);

            grp[S].torsion_dof(2,4,7,8);

            //Threonine
            grp[T].CommonName("Threonine");

            grp[T].init(std::string(" N  , H  , CA , HA , CB , HB , OG1, CG2,") +
                        std::string("1HG2,2HG2,3HG2, HG1, C  , O  "));

            grp[T].set_type("natural aminoacid");

            grp[T].TLC("THR");

            charcode[T]='T';

            olcof["THR"]=T;

            olcof["Threonine"]=T;

            olcof["T"]=T;

            grp[T].torsion_dof(0,2,4,6);

            grp[T].torsion_dof(2,4,6,11);

            grp[T].torsion_dof(2,4,7,8);

            //Cysteine
            grp[C].CommonName("Cysteine");

            grp[C].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , SG ,") +
                        std::string(" HG , C  , O  "));

            grp[C].set_type("natural aminoacid");

            grp[C].TLC("CYS");

            charcode[C]='C';

            olcof["CYS"]=C;

            olcof["Cysteine"]=C;

            olcof["C"]=C;

            grp[C].torsion_dof(0,2,4,7);

            grp[C].torsion_dof(2,4,7,8);

            //Methionine
            grp[M].CommonName("Methionine");

            grp[M].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string("1HG ,2HG , SD , CE ,1HE ,2HE ,3HE , C  ,") +
                        std::string(" O  "));

            grp[M].set_type("natural hydrophobic aminoacid");

            grp[M].TLC("MET");

            charcode[M]='M';

            olcof["MET"]=M;

            olcof["Methionine"]=M;

            olcof["M"]=M;

            grp[M].torsion_dof(0,2,4,7);

            grp[M].torsion_dof(2,4,7,10);

            grp[M].torsion_dof(4,7,10,11);

            grp[M].torsion_dof(7,10,11,12);

            //Proline
            grp[P].CommonName("Proline");

            grp[P].init(std::string(" N  , CA , HA , CB ,1HB ,2HB , CG ,1HG ,") +
                        std::string("2HG , CD ,1HD ,2HD , C  , O  "));

            grp[P].set_type("natural hydrophobic aminoacid");

            grp[P].TLC("PRO");

            charcode[P]='P';

            olcof["PRO"]=P;

            olcof["Proline"]=P;

            olcof["P"]=P;

            //Aspartic_acid
            grp[D].CommonName("Aspartic_acid");

            grp[D].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" OD1, OD2, C  , O  "));

            grp[D].set_type("natural charged aminoacid");

            grp[D].TLC("ASP");

            charcode[D]='D';

            olcof["ASP"]=D;

            olcof["Aspartic_acid"]=D;

            olcof["D"]=D;

            grp[D].torsion_dof(0,2,4,7);

            grp[D].torsion_dof(2,4,7,8);

            //Asparagine
            grp[N].CommonName("Asparagine");

            grp[N].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" OD1, ND2,1HD2,2HD2, C  , O  "));

            grp[N].set_type("natural aminoacid");

            grp[N].TLC("ASN");

            charcode[N]='N';

            olcof["ASN"]=N;

            olcof["Asparagine"]=N;

            olcof["N"]=N;

            grp[N].torsion_dof(0,2,4,7);

            grp[N].torsion_dof(2,4,7,8);

            //Glutamic_acid
            grp[E].CommonName("Glutamic_acid");

            grp[E].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string("1HG ,2HG , CD , OE1, OE2, C  , O  "));

            grp[E].set_type("natural charged aminoacid");

            grp[E].TLC("GLU");

            charcode[E]='E';

            olcof["GLU"]=E;

            olcof["Glutamic_acid"]=E;

            olcof["E"]=E;

            grp[E].torsion_dof(0,2,4,7);

            grp[E].torsion_dof(2,4,7,10);

            grp[E].torsion_dof(4,7,10,11);

            //Glutamine
            grp[Q].CommonName("Glutamine");

            grp[Q].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string("1HG ,2HG , CD , OE1, NE2,1HE2,2HE2, C  ,") +
                        std::string(" O  "));

            grp[Q].set_type("natural aminoacid");

            grp[Q].TLC("GLN");

            charcode[Q]='Q';

            olcof["GLN"]=Q;

            olcof["Glutamine"]=Q;

            olcof["Q"]=Q;

            grp[Q].torsion_dof(0,2,4,7);

            grp[Q].torsion_dof(2,4,7,10);

            grp[Q].torsion_dof(4,7,10,11);

            //Lysine
            grp[K].CommonName("Lysine");

            grp[K].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string("1HG ,2HG , CD ,1HD ,2HD , CE ,1HE ,2HE ,") +
                        std::string(" NZ ,1HZ ,2HZ ,3HZ , C  , O  "));

            grp[K].set_type("natural charged aminoacid");

            grp[K].TLC("LYS");

            charcode[K]='K';

            olcof["LYS"]=K;

            olcof["Lysine"]=K;

            olcof["K"]=K;

            grp[K].torsion_dof(0,2,4,7);

            grp[K].torsion_dof(2,4,7,10);

            grp[K].torsion_dof(4,7,10,13);

            grp[K].torsion_dof(7,10,13,16);

            grp[K].torsion_dof(10,13,16,17);

            //Arginine
            grp[R].CommonName("Arginine");

            grp[R].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string("1HG ,2HG , CD ,1HD ,2HD , NE , HE , CZ ,") +
                        std::string(" NH1,1HH1,2HH1, NH2,1HH2,2HH2, C  , O  "));

            grp[R].set_type("natural charged aminoacid");

            grp[R].TLC("ARG");

            charcode[R]='R';

            olcof["ARG"]=R;

            olcof["Arginine"]=R;

            olcof["R"]=R;

            grp[R].torsion_dof(0,2,4,7);

            grp[R].torsion_dof(2,4,7,10);

            grp[R].torsion_dof(4,7,10,13);

            grp[R].torsion_dof(7,10,13,15);

            //Histidine
            grp[H].CommonName("Histidine");

            grp[H].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" CD2, NE2, CE1, ND1, HD2, HE1, HD1, C  ,") +
                        std::string(" O  "));

            grp[H].set_type("natural aminoacid");

            grp[H].TLC("HIS");

            charcode[H]='H';

            olcof["HIS"]=H;

            olcof["Histidine"]=H;

            olcof["H"]=H;

            grp[H].torsion_dof(0,2,4,7);

            grp[H].torsion_dof(2,4,7,8);

            //Phenylalanine
            grp[F].CommonName("Phenylalanine");

            grp[F].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" CD1, CE1, CZ , CE2, CD2, HD1, HE1, HZ ,") +
                        std::string(" HE2, HD2, C  , O  "));

            grp[F].set_type("natural hydrophobic aminoacid");

            grp[F].TLC("PHE");

            charcode[F]='F';

            olcof["PHE"]=F;

            olcof["Phenylalanine"]=F;

            olcof["F"]=F;

            grp[F].torsion_dof(0,2,4,7);

            grp[F].torsion_dof(2,4,7,8);

            //Tyrosine
            grp[Y].CommonName("Tyrosine");

            grp[Y].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" CD1, CE1, CZ , CE2, CD2, HD1, HE1, OH ,") +
                        std::string(" HE2, HD2, HH , C  , O  "));

            grp[Y].set_type("natural hydrophobic aminoacid");

            grp[Y].TLC("TYR");

            charcode[Y]='Y';

            olcof["TYR"]=Y;

            olcof["Tyrosine"]=Y;

            olcof["Y"]=Y;

            grp[Y].torsion_dof(0,2,4,7);

            grp[Y].torsion_dof(2,4,7,8);

            grp[Y].torsion_dof(9,10,15,18);

            //Tryptophan
            grp[W].CommonName("Tryptophan");

            grp[W].init(std::string(" N  , H  , CA , HA , CB ,1HB ,2HB , CG ,") +
                        std::string(" CD1, NE1, CE2, CD2, CE3, CZ3, CH2, CZ2,") +
                        std::string(" HD1, HE1, HE3, HZ3, HH2, HZ2, C  , O  "));

            grp[W].set_type("natural hydrophobic aminoacid");

            grp[W].TLC("TRP");

            charcode[W]='W';

            olcof["TRP"]=W;

            olcof["Tryptophan"]=W;

            olcof["W"]=W;

            grp[W].torsion_dof(0,2,4,7);

            grp[W].torsion_dof(2,4,7,8);

            //D-Proline
            grp[DPR].CommonName("D-Proline");

            grp[DPR].init(std::string(" N  , CA , HA , CB ,1HB ,2HB , CG ,1HG ,") +
                          std::string("2HG , CD ,1HD ,2HD , C  , O  "));

            grp[DPR].set_type("synthetic hydrophobic aminoacid");

            grp[DPR].TLC("DPR");

            charcode[DPR]='p';

            olcof["DPR"]=DPR;

            olcof["D-Proline"]=DPR;

            //Acetyl
            grp[ACE].CommonName("Acetyl");

            grp[ACE].init(std::string(" C  , O  , CH3,1H  ,2H  ,3H  "));

            grp[ACE].set_type("N-terminal endgroup");

            grp[ACE].TLC("ACE");

            charcode[ACE]='X';

            olcof["ACE"]=ACE;

            olcof["Acetyl"]=ACE;

            grp[ACE].torsion_dof(6,0,2,3);

            //Amide
            grp[NH2].CommonName("Amide");

            grp[NH2].init(std::string(" N  ,1HN ,2HN "));

            grp[NH2].set_type("C-terminal endgroup");

            grp[NH2].TLC("NH2");

            charcode[NH2]='X';

            olcof["NH2"]=NH2;

            olcof["Amide"]=NH2;

            //Succinyl
            grp[SUC].CommonName("Succinyl");

            grp[SUC].init(std::string(" CO , O  , C3 ,1H3 ,2H3 , C2 ,1H2 ,2H2 ,") +
                          std::string(" CO2, O11, O12"));

            grp[SUC].set_type("N-terminal endgroup");

            grp[SUC].TLC("SUC");

            charcode[SUC]='X';

            olcof["SUC"]=SUC;

            olcof["Succinyl"]=SUC;

            grp[SUC].torsion_dof(11,0,2,3);

            grp[SUC].torsion_dof(0,2,5,6);

            grp[SUC].torsion_dof(2,5,8,9);

            //NMethyl
            grp[NME].CommonName("NMethyl");

            grp[NME].init(std::string(" N  , H  , CH3,1H  ,2H  ,3H  "));

            grp[NME].set_type("C-terminal endgroup");

            grp[NME].TLC("NME");

            charcode[NME]='X';

            olcof["NME"]=NME;

            olcof["NMethyl"]=NME;

            grp[NME].torsion_dof(-1,0,2,3);

            //VoidEG
            grp[VOIDEG].CommonName("VoidEG");

            grp[VOIDEG].init(std::string(""));

            grp[VOIDEG].set_type("Blank endgroup");

            grp[VOIDEG].TLC("___");

            charcode[VOIDEG]='X';

            olcof["___"]=VOIDEG;

            olcof["VoidEG"]=VOIDEG;
            already_initialized=true;
        }
    }
}

