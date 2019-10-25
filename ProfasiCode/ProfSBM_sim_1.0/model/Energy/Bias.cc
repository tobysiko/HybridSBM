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

#include "Bias.hh"
#include <cmath>

using namespace std;
using UnivConstants::pi;

namespace prf
{
    Bias::Bias() : Energy()
    {
        Name("Bias");
        //default parameters
        kappaH=1.2; kappaO=1.2; kcou=6.0;
        qC=.42; qO=-.42; qN=-.20; qH=.20;
        grdtyp=2;
    }

    Bias::~Bias() {}

    void Bias::init()
    {
        if (initialized) return;
        Mebias.resize(NC());
        Vebias.resize(NC());

        //The array qPos[][][] will list identities and q[][][]
        //charges of all peptide unit atoms. q[i][j][] should contain
        //Cj-1, Oj-1, Nj and Hj, in that order. If there are no
        //capping groups, bias will not be calculated for the first
        //and last residues, e.g. q[i][0][] contains the c-terminus
        //peptide group if there is any, otherwise it is empty.
        q.resize(NC());
        qPos.resize(NC());

        for (int i=0; i<NC(); ++i) {
            q[i].resize(p->Chain(i)->numAminoAcids()+1);
            qPos[i].resize(q[i].size());
            Mebias[i].resize(p->Chain(i)->numAminoAcids(),0);
            Vebias[i].resize(p->Chain(i)->numAminoAcids(),0);

            //handle inner groups first

            for (int j=1; j<p->Chain(i)->numAminoAcids(); ++j) {
                AminoAcid* prev = p->Chain(i)->AA(j-1);
                AminoAcid* curr = p->Chain(i)->AA(j);

                //assign C and 0 of peptide unit j
                qPos[i][j].push_back(prev->Cprime().UniqueId());
                q[i][j].push_back(qC);
                qPos[i][j].push_back(prev->Oc().UniqueId());
                q[i][j].push_back(qO);

                //assign N and H of peptide unit j (if applicable)

                if (curr->OLC()!=P && curr->OLC()!=DPR) {
                    qPos[i][j].push_back(curr->Nitrogen().UniqueId());
                    q[i][j].push_back(qN);
                    qPos[i][j].push_back(curr->Nitrogen().UniqueId()+1);
                    q[i][j].push_back(qH);
                }
            }

            //handle N-terminus
            EndGroup* eg;

            if ((eg=p->Chain(i)->NtermLigand()) != NULL) {
                if (eg->pep_bond_link()) {
                    int ic=0,io=1;
                    eg->pep_bond_atoms(ic,io);
                    qPos[i][0].push_back(ic);
                    q[i][0].push_back(qC);
                    qPos[i][0].push_back(io);
                    q[i][0].push_back(qO);
                }

                AminoAcid* ntres = p->Chain(i)->AA(0);

                if (ntres->OLC()!=P && ntres->OLC()!=DPR) {
                    qPos[i][0].push_back(ntres->Nitrogen().UniqueId());
                    q[i][0].push_back(qN);
                    qPos[i][0].push_back(ntres->Nitrogen().UniqueId()+1);
                    q[i][0].push_back(qH);
                }
            }

            if ((eg=p->Chain(i)->CtermLigand()) != NULL) {
                int jlast=p->Chain(i)->numAminoAcids()-1;
                AminoAcid* ctres = p->Chain(i)->AA(jlast);
                qPos[i][jlast+1].push_back(ctres->Cprime().UniqueId());
                q[i][jlast+1].push_back(qC);
                qPos[i][jlast+1].push_back(ctres->Cprime().UniqueId()+1);
                q[i][jlast+1].push_back(qO);

                if (eg->pep_bond_link()) {
                    int in=0,ih=1;
                    eg->pep_bond_atoms(in,ih);
                    qPos[i][jlast+1].push_back(in);
                    q[i][jlast+1].push_back(qN);
                    qPos[i][jlast+1].push_back(ih);
                    q[i][jlast+1].push_back(qH);
                }
            }
        }

        changed.resize(p->NumberOfLigands());

        nchanges=0;
        initialized=true;
    }

    double Bias::evaluate()
    {
        vval=0;

        for (int ich=0;ich<NC();++ich) {
            int N=p->Chain(ich)->numAminoAcids();

            for (int i=0;i<N;++i) {
                double eterm=biasTerm(ich,i);
                vval+=(Mebias[ich][i]=Vebias[ich][i]=eterm);
            }
        }

        delv=0;

        return vval;
    }

    double Bias::gradientXYZ(std::valarray<double> &gx)
    {
        vval=0;

        for (int ich=0;ich<NC();++ich) {
            int N=p->Chain(ich)->numAminoAcids();

            for (int i=0;i<N;++i) {
                double eterm=coulomb_with_grd(ich,i,gx)+
                             repulsion_with_grd(ich,i,gx);
                vval+=eterm;
            }
        }

        delv=0;

        return vval;
    }

    double Bias::deltaE(Update *updt)
    {
        delv=0;
        nchanges=0;

        if (updt->sidechain_update()||
            updt->rigid_chain_update()) return delv=0;

        double eterm=0,de=0;

        int numsets=updt->n_residue_rigid_ranges();

        for (int i=0;i<numsets;++i) {
            int i1=0,i2=0;
            updt->residue_rigid_range(i,i1,i2);

            if ((i2-i1)==0) {
                Ligand *lg=p->ligand(i1);

                if (lg->isAA()) {
                    i2=lg->LocatedOn();
                    i1=lg->UniqueId()-p->Chain(i2)->AA(0)->UniqueId();
                    eterm=biasTerm(i2,i1);

                    if ((de=eterm-Mebias[i2][i1])!=0) {
                        changed[nchanges].first=i2;
                        changed[nchanges++].second=i1;
                        Vebias[i2][i1]=eterm;
                    }

                    delv+=de;
                }
            }
        }

        return delv;
    }

    void Bias::Accept(Update *updt)
    {
        vval+=delv;

        for (; nchanges; --nchanges)
            Mebias[changed[nchanges-1].first][changed[nchanges-1].second]=
                Vebias[changed[nchanges-1].first][changed[nchanges-1].second];
    }

    void Bias::rangeEstimate(double &x1, double &x2)
    {
        int npro=0,ngly=0,nAA=p->NumberOfResidues(),nendgrp=0;

        for (int j=0;j<NC();++j) {
            for (int i=0;i<p->Chain(j)->numAminoAcids();++i) {
                prf::OneLetterCode olc=p->Chain(j)->AA(i)->OLC();
                if (olc==P || olc==DPR) ++npro;
                else if (olc==G) ++ngly;
            }

            if (p->Chain(j)->NtermLigand()!=NULL) ++nendgrp;

            if (p->Chain(j)->CtermLigand()!=NULL) ++nendgrp;
        }

        x1=-0.6*(nAA-npro+nendgrp-2*NC());

        x2=0.25*(nAA-ngly-npro)*kappaH+0.25*(nAA-ngly)*kappaO;
        //There is no good reason behind these ranges. They are just
        //preliminary estimates based on experience.
    }

    double Bias::biasTerm(int ich, int ires)
    {
        if (ires < 0) return 0;

//         return repulsion(ich,ires)+coulomb(ich,ires);
        return coulomb(ich,ires)+repulsion(ich,ires);
    }

    double Bias::repulsion(int ich, int ires)
    {
        double e=0;
        //if there is only one peptide group in connection with the
        //residue in question the energy will be zero

        if (qPos[ich][ires].size()<2 || qPos[ich][ires+1].size()<2)
            return 0;

        //skip Glycine
        if (p->Chain(ich)->AA(ires)->OLC()==G)
            return 0;

        //H H repulsion (if no Proline involved)
        if (qPos[ich][ires].size()==4 && qPos[ich][ires+1].size()==4) {
            int in=qPos[ich][ires][2], ih=qPos[ich][ires][3];
            int inn=qPos[ich][ires+1][2], ihn=qPos[ich][ires+1][3];
            e+=hrf(ih,inn,ihn,in,kappaH);
        }

        //O O repsulsion
        int icp=qPos[ich][ires][0], iop=qPos[ich][ires][1];
        int ic=qPos[ich][ires+1][0], io=qPos[ich][ires+1][1];

        e+=hrf(iop,ic,io,icp,kappaO);

        return e;
    }

    double Bias::repulsion_with_grd(int ich, int ires,
                                    std::valarray<double> &gx)
    {
        double e=0;
        //if there is only one peptide group in connection with the
        //residue in question the energy will be zero

        if (qPos[ich][ires].size()<2 || qPos[ich][ires+1].size()<2)
            return 0;

        //skip Glycine
        if (p->Chain(ich)->AA(ires)->OLC()==G)
            return 0;

        //H H repulsion (if no Proline involved)
        if (qPos[ich][ires].size()==4 && qPos[ich][ires+1].size()==4) {
            int in=qPos[ich][ires][2], ih=qPos[ich][ires][3];
            int inn=qPos[ich][ires+1][2], ihn=qPos[ich][ires+1][3];
            e+=dhrf(ih,inn,ihn,in,kappaH,gx);
        }

        //O O repsulsion
        int icp=qPos[ich][ires][0], iop=qPos[ich][ires][1];
        int ic=qPos[ich][ires+1][0], io=qPos[ich][ires+1][1];

        e+=dhrf(iop,ic,io,icp,kappaO,gx);

        return e;
    }

    double Bias::hrf(int i1, int i2, int i3, int i4, double strk)
    {
        double d=0;
        if ((d=min(AtomCoordinates::s(i1,i2),
                   AtomCoordinates::s(i3,i4))-
             AtomCoordinates::s(i1,i3)) > 0)
            return strk*tanh(3*d);
        else return 0;
    }

    double Bias::dhrf(int i1, int i2, int i3, int i4, double strk,
                      std::valarray<double> &gx)
    {
        Vector3 r12=AtomCoordinates::vec(i1)-AtomCoordinates::vec(i2);
        Vector3 r34=AtomCoordinates::vec(i3)-AtomCoordinates::vec(i4);
        Vector3 r13=AtomCoordinates::vec(i1)-AtomCoordinates::vec(i3);

        double s12=r12.mag(),s34=r34.mag(),s13=r13.mag();
        double d=min(s12,s34)-s13;
        if (d<=0) return 0;
        double grdmg=cosh(3*d);
        grdmg=3*strk/(grdmg*grdmg);
        r13*=(-grdmg/s13);
        gx[3*i1]+=r13.x();
        gx[3*i1+1]+=r13.y();
        gx[3*i1+2]+=r13.z();
        gx[3*i3]-=r13.x();
        gx[3*i3+1]-=r13.y();
        gx[3*i3+2]-=r13.z();
        if (s12<s34) {
            r12*=(grdmg/s12);
            gx[3*i1]+=r12.x();
            gx[3*i1+1]+=r12.y();
            gx[3*i1+2]+=r12.z();
            gx[3*i2]-=r12.x();
            gx[3*i2+1]-=r12.y();
            gx[3*i2+2]-=r12.z();
        } else {
            r34*=(grdmg/s34);
            gx[3*i3]+=r34.x();
            gx[3*i3+1]+=r34.y();
            gx[3*i3+2]+=r34.z();
            gx[3*i4]-=r34.x();
            gx[3*i4+1]-=r34.y();
            gx[3*i4+2]-=r34.z();
        }
        return strk*tanh(3*d);
    }

    double Bias::coulomb(int ich,int ires)
    {
        double e=0;
        //if there are no end-groups q[ich][0] and q[ich][NAA] are
        //empty, i.e. the loop is never entered
        //Resiude ires is between peptide bonds ires and ires+1.

        for (size_t i=0;i<q[ich][ires].size();++i) {
            for (size_t j=0;j<q[ich][ires+1].size();++j) {
                e+=q[ich][ires][i]*q[ich][ires+1][j]/
                   AtomCoordinates::s(qPos[ich][ires][i],
                                      qPos[ich][ires+1][j]);
            }
        }

        return kcou*e;
    }

    double Bias::coulomb_with_grd(int ich,int ires,
                                  std::valarray<double> &gx)
    {
        double e=0,r=1,et=0;
        Vector3 v12;
        int iat=0,jat=0;
        //if there are no end-groups q[ich][0] and q[ich][NAA] are
        //empty, i.e. the loop is never entered
        //Resiude ires is between peptide bonds ires and ires+1.

        for (size_t i=0;i<q[ich][ires].size();++i) {
            for (size_t j=0;j<q[ich][ires+1].size();++j) {
                v12=AtomCoordinates::vec(iat=qPos[ich][ires][i])
                    -AtomCoordinates::vec(jat=qPos[ich][ires+1][j]);
                r=v12.mag();
                et=q[ich][ires][i]*q[ich][ires+1][j]/r;
                e+=et;
                v12*=(-et/r/r);
                gx[3*iat]+=v12.x();
                gx[3*iat+1]+=v12.y();
                gx[3*iat+2]+=v12.z();
                gx[3*jat]-=v12.x();
                gx[3*jat+1]-=v12.y();
                gx[3*jat+2]-=v12.z();
            }
        }

        return kcou*e;
    }
}

