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

#include "ChargedSCInteraction.hh"

using namespace std;

namespace prf
{
    const int ChargedSCInteraction::nCHgrp=7;
    const double ChargedSCInteraction::hpstr[]={ -1, -1, 1, 1, -1, 1, -1 };
    const int cdnt=5;
    const int cdct=6;
    string ChargedSCInteraction::hatms[][4]={
        {" OD1"," OD2","    ","    "}, //D
        {" OE1"," OE2","    ","    "}, //E
        {"1HZ ","2HZ ","3HZ ","    "}, //K
        {" NE "," CZ "," NH1"," NH2"}, //R
        {" O11"," O12","    ","    "}, //SUC
        {"1H  ","2H  ","3H  ","    "}, //N-terminal end of a chain
        {" O  "," OXT","    ","    "}  //C-terminal end of a chain
    };

    int ChargedSCInteraction::htyp(OneLetterCode cod)
    {
        OneLetterCode hpcdord[]={NONE,D,E,K,R,SUC};
        int i=0;

        for (i=5;i>0&&hpcdord[i]!=cod;--i);

        return i-1;
    }

    ChargedSCInteraction::ChargedSCInteraction() : Energy()
    {
        Name("ChargedSCInteraction");
        a2=3.7*3.7;
        b2=4.5*4.5;
        d2=b2-a2;
        slpe=-1.0/d2;
        strnn=1.0;
        strnnn=1.0;
        sclfct=1.0;
        grdtyp=2;
    }

    ChargedSCInteraction::~ChargedSCInteraction() {}

    double ChargedSCInteraction::InterChain(int ich, int jch)
    {
        double e=0;
        int hi,hj,jbeg;

        for (int i=p->chain_start(ich);i<p->chain_end(ich);++i) {
            if ((hi=ih[i])<0) continue;

            if (ich==jch) jbeg=i+1;
            else jbeg=p->chain_start(jch);

            for (int j=jbeg;j<p->chain_end(jch);++j) {
                if ((hj=ih[j])<0) continue;

                e+=Mehp[NCAA*hi+hj];
            }
        }

        return sclfct*e;
    }

    double ChargedSCInteraction::intercg()
    {
        double e=0;

        for (int ich=0;ich<NC();++ich) {
            for (int jch=ich+1;jch<NC();++jch) {
                e+=InterChain(ich,jch);
            }
        }

        return e;
    }

    void ChargedSCInteraction::init()
    {
        if (initialized) return;
        Logger blog;
        double tmpx;
        N=p->NumberOfLigands(),NCAA=0;

        for (int ich=0;ich<p->NumberOfChains();++ich) {
            if (p->Chain(ich)->AA(0)->hasNTerminal()) ++NCAA;
            if (p->Chain(ich)->AA(p->Chain(ich)->numAminoAcids()-1)->hasCTerminal()) ++NCAA;
        }

        for (int ia=0;ia<p->NumberOfLigands();++ia) {
            if (htyp(p->ligand(ia)->OLC())>=0) ++NCAA;
        }

        blog(15)<<Name()<<"> Counted "<<NCAA<<" charged side chains\n";

        //NCAA = total number of  charged amino acids or endgroups
        //In the following, the strange initialisations are an attempt to
        //try and crash the program if something illegal happens.
        strength.resize(NCAA*NCAA,0.0);
        atom.resize(NCAA*4,-200000);
        natom.resize(NCAA,0);
        std::vector<int> tmpst(NCAA,0);
        jh.resize(NCAA,-1);
        ih.resize(N,-1);
        nh.resize(N,0);
        seq.resize(N,NONE);
        chainof.resize(NCAA,-1);
        Mehp.resize(NCAA*NCAA,0);
        Vehp.resize(NCAA*NCAA,0);
        changed.resize(NCAA*NCAA,0);
        Ligand *lg;


        for (int i1=0,k1=0,l1=0;i1<NC();++i1) {
            AminoAcid *aa0=p->Chain(i1)->AA(0);

            if (aa0->hasNTerminal()) {
                ++nh[k1];
                chainof[l1]=i1;

                for (int i=0;i<4;++i) {
                    Atom *at1=aa0->labeled_atom(hatms[cdnt][i]);

                    if (at1!=NULL) {atom[4*l1+i]=at1->UniqueId();natom[l1]++;}
                }

                tmpst[l1]=hpstr[cdnt];

                if (ih[k1]==-1) ih[k1]=l1;

                jh[l1]=k1;

                l1++;
            }

            for (int j1=0;j1<p->Chain(i1)->numLigands();++j1) {
                lg=p->Chain(i1)->memberLigand(j1);

                if (htyp(seq[k1]=lg->OLC()) !=-1) {
                    chainof[l1]=i1;
                    ++nh[k1];

                    for (int i=0;i<4;++i) {
                        Atom *at1=lg->labeled_atom(hatms[htyp(seq[k1])][i]);

                        if (at1!=NULL) {atom[4*l1+i]=at1->UniqueId();natom[l1]++;}
                    }

                    tmpst[l1]=hpstr[htyp(seq[k1])];

                    if (ih[k1]==-1) ih[k1]=l1;

                    jh[l1]=k1;

                    l1++;
                }

                if (j1<(p->Chain(i1)->numLigands()-1)) ++k1;
            }

            AminoAcid *aa1=p->Chain(i1)->AA(p->Chain(i1)->numAminoAcids()-1);

            if (aa1->hasCTerminal()) {
                ++nh[k1];
                chainof[l1]=i1;

                for (int i=0;i<4;++i) {
                    Atom *at1=aa1->labeled_atom(hatms[cdct][i]);

                    if (at1!=NULL) {atom[4*l1+i]=at1->UniqueId();natom[l1]++;}
                }

                tmpst[l1]=hpstr[cdct];

                if (ih[k1]==-1) ih[k1]=l1;

                jh[l1]=k1;

                l1++;
            }
            ++k1;
        }

        for (int l1=0;l1<NCAA;++l1) {
            for (int l2=0;l2<NCAA;++l2) {
                tmpx=1.5*tmpst[l1]*tmpst[l2];

                if (chainof[l1]==chainof[l2]) {

                    if (jh[l2]==(jh[l1]+1)||jh[l2]==(jh[l1]-1)) tmpx*=strnn;

                    if (jh[l2]==(jh[l1]+2)||jh[l2]==(jh[l1]-2)) tmpx*=strnnn;
                }

                strength[l1*NCAA+l2]=tmpx;
            }
        }

        blog(50)<<"List of Ligands participating in charged side chain interactions ...\n";

        for (int i=0;i<NCAA;++i) {
            blog<<i<<"  "<<Groups::mapOLC2Char(seq[jh[i]])<<" on chain "
            <<chainof[i]<<" index "<<p->ligand(jh[i])->SeqSerial()
            <<"\n";
        }

        blog<<"Pairwise strength matrix for charged side chain interactions ...\n";

        for (int i=0;i<NCAA;++i) {
            for (int j=0;j<NCAA;++j) {
                blog<<strength[i*NCAA+j]<<"  ";
            }

            blog<<"\n";
        }


        double max_at_dist=(14+sqrt(a2));

        max_at_dist2=max_at_dist*max_at_dist;
        max_at_dist=AtomCoordinates::boxL()-max_at_dist;
        max_cmp_dst2=max_at_dist*max_at_dist;
        initialized=true;
    }

    void ChargedSCInteraction::DisplayMatrix()
    {
        for (int i=0;i<NCAA;++i) {
            for (int j=0;j<NCAA;++j) prf::cout <<Mehp[NCAA*i+j]<<"\t";

            prf::cout <<"\n";
        }

        for (int i=0;i<NCAA;++i) {
            for (int j=0;j<NCAA;++j) prf::cout <<Vehp[NCAA*i+j]<<"\t";

            prf::cout <<"\n";
        }
    }

    double ChargedSCInteraction::evaluate()
    {
        double eterm=0;
        vval=delv=0;
        nchanges=0;

        for (int i=0;i<NCAA;++i) {
            for (int j=i+1;j<NCAA;++j) {
                eterm=cg_pair(i,j,chainof[i]==chainof[j]?
                              AtomCoordinates::s2:
                              AtomCoordinates::dist2);
                vval+=eterm;
                /*                if (fabs(eterm)>1e-6) {
                                    prf::cout<<"contribution "<<i<<", "<<j<<" = "<<eterm<<" with "<<strength[i*NCAA+j]<<"\n";
                                }*/
                Mehp[NCAA*i+j]=Mehp[NCAA*j+i]=eterm;
            }
        }

        return vval*=sclfct;
    }

    double ChargedSCInteraction::gradientXYZ(std::valarray<double> &gx)
    {
        double eterm=0;
        vval=delv=0;
        nchanges=0;
        gx=0;
        for (int i=0;i<NCAA;++i) {
            for (int j=i+1;j<NCAA;++j) {
                eterm=dcg_pair(i,j,chainof[i]==chainof[j],gx);
                vval+=eterm;
                Mehp[NCAA*i+j]=Mehp[NCAA*j+i]=eterm;
            }
        }

        return vval*=sclfct;
    }

    void ChargedSCInteraction::rangeEstimate(double &x1, double &x2)
    {
        x1=-1.0*sclfct*max(NCAA,1);

        x2=-0.01*x1;
        //mere guess.
    }

    double ChargedSCInteraction::deltaE(Update *updt)
    {
        delv=0;
        nchanges=0;

        if (updt->sidechain_update()) {
            int ind=ih[updt->change(0).info.group];

            if (ind<0) return 0;
        }

        double eterm=0,de=0;

        updt->n_residue_rigid_ranges();
        vector<pair<int,int> > *rng=updt->residue_rigid_ranges();

        for (int it=0; it<updt->n_residue_rigid_ranges();++it) {
            int ires1=(*rng)[it].first,   ires2=(*rng)[it].second+1;
//             prf::cout<<"Update "<<updt->Name()<<" res range "<<ires1<<" to "<<ires2<<"\n";
            int ihres=0,jhres=0;

            for (int ires=ires1;ires<ires2;++ires) {
                if ((ihres=ih[ires])<0) continue;

                for (;ihres<ih[ires]+nh[ires];++ihres) {
                    for (int jres=0,prng=0;jres<N;) {
                        if ((jhres=ih[jres])<0||(jres>=ires1 && jres<ires2 && jres<ires)) {++jres; continue;}

                        if (prng<it && jres>=(*rng)[prng].first) {
                            jres=(*rng)[prng++].second+1;
                            continue;
                        }

                        for (; jhres < ih[jres]+nh[jres];++jhres) {
                            if (jhres==ihres || (jres==ires && jhres<ihres)) continue;

                            if (chainof[ihres]==chainof[jhres])
                                eterm=cg_pair(ihres,jhres,AtomCoordinates::s2);
                            else eterm=cg_pair(ihres,jhres,AtomCoordinates::dist2);

                            delv+=(de=eterm-Mehp[ihres*NCAA+jhres]);

//         prf::cout<<"new value for "<<ihres<<", "<<jhres<<", "<<ires<<", "<<jres<<" = "<<eterm<<", "<<de<<"\n";
                            if (de!=0) {
                                Vehp[changed[nchanges++]=ihres*NCAA+jhres]=eterm;
                                Vehp[changed[nchanges++]=jhres*NCAA+ihres]=eterm;
                            }
                        }

                        ++jres;
                    }
                }
            }
        }

        return (delv*=sclfct);
    }

    void ChargedSCInteraction::Accept(Update *updt)
    {
        for (;nchanges;--nchanges)
            Mehp[changed[nchanges-1]]=Vehp[changed[nchanges-1]];

        vval+=delv;
    }

    double ChargedSCInteraction::cg_pair(int i, int j, double(*distf)(int,int))
    {
        /* i and j are indices for NCAA sized arrays. */
        double ans=strength[i*NCAA+j]*contact_frac(i,j,distf);
        return ans;
    }

    double ChargedSCInteraction::cg_contact_frac(int iaa, int jaa,
            double(*distf)(int,int))
    {
        /* iaa and jaa are amino acid indices relative to the
        entire system. This function is intended only for use
        in connection with observables.*/
        int i=ih[iaa],j=ih[jaa];

        if (i<0||j<0) return 0.0;

        return contact_frac(i,j,distf);
    }

    double ChargedSCInteraction::contact_frac(int i, int j,
            double(*distf)(int,int))
    {
        /* i and j are indices for NCAA sized arrays. */
        int ni,nj,k,l,m,n;
        double tmp,sum=0;
        ni=natom[i];
        nj=natom[j];

        for (n=0;n<ni+nj;n++) r2min[n]=1e10;

        for (n=0;n<ni;n++) {
            k=atom[i*4+n];

            for (m=0;m<nj;m++) {
                l=atom[j*4+m];
                tmp=distf(k,l);

                if (tmp>max_at_dist2 && tmp<max_cmp_dst2) return 0;

                if (tmp<r2min[n]) r2min[n]=tmp;

                if (tmp<r2min[ni+m]) r2min[ni+m]=tmp;
            }
        }

        for (n=0;n<ni+nj;n++) {
            if (r2min[n]>b2) continue;

            if (r2min[n]<a2) sum++;
            else sum+= (b2-r2min[n]) /d2;
        }

        /*                if (fabs(sum)>1e-6) {
                            prf::cout<<"contribution "<<i<<", "<<j<<" = "<<sum<<"\n";
              prf::cout<<"atoms for group "<<i<<": \n";
              for (int ia=0;ia<natom[i];++ia) {
           prf::cout<<atom[4*i+ia]<<"  ";
              }
              prf::cout<<"\n";
              prf::cout<<"atoms for group "<<j<<": \n";
              for (int ia=0;ia<natom[j];++ia) {
           prf::cout<<atom[4*j+ia]<<"  ";
              }
              prf::cout<<"\n";
                        }*/
        return std::min(1.0,sum/(ni+nj));
    }
}

double ChargedSCInteraction::dcg_pair(int i, int j, bool perdhint,
                                      std::valarray<double> &gx)
{
    double pstr=0;
    int ni,nj,k,l,m,n;
    double tmp,sum=0;
    ni=natom[i];
    nj=natom[j];
    pstr=strength[i*NCAA+j]*sclfct;
    for (n=0;n<ni+nj;n++) r2min[n]=1e10;

    for (n=0;n<ni;n++) {
        k=atom[i*4+n];

        for (m=0;m<nj;m++) {
            l=atom[j*4+m];
            tmp=perdhint?AtomCoordinates::dist2(k,l):AtomCoordinates::s2(k,l);

            if (tmp>max_at_dist2 && tmp<max_cmp_dst2) return 0;

            if (tmp<r2min[n]) {
                r2min[n]=tmp;
                r2minpartner[n]=m;
            }

            if (tmp<r2min[ni+m]) {
                r2min[ni+m]=tmp;
                r2minpartner[ni+m]=n;
            }
        }
    }

    for (n=0;n<ni+nj;++n) {
        if (r2min[n]>b2) continue;
        if (r2min[n]<a2) sum++;
        else sum+=(b2-r2min[n])/d2;
    }

    if (sum<=0) return 0;
    if (sum<(ni+nj)) {
        pstr/=(ni+nj);
        for (n=0;n<ni+nj;++n) {
            if (r2min[n]>b2 or r2min[n]<a2) continue;
            if (n<ni) {
                m=atom[i*4+n];
                l=atom[4*j+r2minpartner[n]];
            } else {
                m=atom[4*j+n-ni];
                l=atom[4*i+r2minpartner[n]];
            }
            Vector3 dirn;
            if (perdhint) dirn=AtomCoordinates::sep(m,l);
            else dirn=AtomCoordinates::diff(m,l);
            dirn*=(2*pstr*slpe);
            gx[3*m]+=dirn.x();
            gx[3*m+1]+=dirn.y();
            gx[3*m+2]+=dirn.z();
            gx[3*l]-=dirn.x();
            gx[3*l+1]-=dirn.y();
            gx[3*l+2]-=dirn.z();
        }
        return pstr*sum;
    } else return pstr;
}
