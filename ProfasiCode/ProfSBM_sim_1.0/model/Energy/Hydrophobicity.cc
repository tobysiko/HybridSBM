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

#include "Hydrophobicity.hh"

using namespace std;

namespace prf
{
    const int Hydrophobicity::nHPgrp=11;
    const double Hydrophobicity::hpstr[]={ 0.3, 0.4, 0.4, 0.6, 0.8, 0.8, 0.8, 0.8, 1.1, 1.6, 1.6 };

    string Hydrophobicity::hatms[][6]={
        {" CB "," CG ","    ","    ","    ","    "}, //R
        {" CB "," CG "," SD "," CE ",    "",    ""}, //M
        {" CB "," CG "," CD ","    ","    ","    "}, //K
        {" CB "," CG1"," CG2",    "",    "",    ""}, //V
        {" CB "," CG1"," CG2"," CD1",    "",    ""}, //I
        {" CB "," CG "," CD1"," CD2",    "",    ""}, //L
        {" CB "," CG "," CD ",    "",    "",    ""}, //P
        {" CB "," CG "," CD ",    "",    "",    ""}, //DPR
        {" CG "," CD1"," CE1"," CZ "," CE2"," CD2"}, //Y
        {" CG "," CD1"," CE1"," CZ "," CE2"," CD2"}, //F
        {" CG "," CD1"," CD2"," CE3"," CZ3"," CH2"}  //W
    };

    int Hydrophobicity::htyp(OneLetterCode cod)
    {
        OneLetterCode hpcdord[]={NONE,R,M,K,V,I,L,P,DPR,Y,F,W};
        int i=0;

        for (i=nHPgrp;i>0&&hpcdord[i]!=cod;--i);

        return i-1;
    }

    Hydrophobicity::Hydrophobicity() : Energy()
    {
        Name("Hydrophobicity");
        a2=3.7*3.7;
        b2=4.5*4.5;
        d2=b2-a2;
        slpe=-1.0/d2;
        strnn=0.0;
        strnnn=0.5;
        sclfct=1.0;
        grdtyp=2;
        focalRes = -1;
        focalRes2 = -1;
        focalResScale = 1.0;
    }

    Hydrophobicity::~Hydrophobicity() {}

    double Hydrophobicity::InterChain(int ich, int jch)
    {
        double e=0;
        int hi,hj,jbeg;

        for (int i=p->chain_start(ich);i<p->chain_end(ich);++i) {
            if ((hi=ih[i])<0) continue;

            if (ich==jch) jbeg=i+1;
            else jbeg=p->chain_start(jch);

            for (int j=jbeg;j<p->chain_end(jch);++j) {
                if ((hj=ih[j])<0) continue;

                e+=Mehp[NHAA*hi+hj];
            }
        }

        return sclfct*e;
    }

    double Hydrophobicity::interhp()
    {
        double e=0;

        for (int ich=0;ich<NC();++ich) {
            for (int jch=ich+1;jch<NC();++jch) {
                e+=InterChain(ich,jch);
            }
        }

        return e;
    }

    void Hydrophobicity::init()
    {
        if (initialized) return;
        Logger blog;
        double tmpx;
        N=p->NumberOfLigands(),NHAA=0;

        for (int i=0;i<p->NumberOfLigands();++i) {
            if (htyp(p->ligand(i)->OLC())>=0) ++NHAA;
        }

        blog(15)<<Name()<<"> Counted "<<NHAA<<" hydrophobic ligands\n";

        //NHAA = total number of  hydrophobic amino acids or endgroups
        //In the following, the strange initialisations are an attempt to
        //try and crash the program if something illegal happens.
        strength.resize(NHAA*NHAA,0.0);
        atom.resize(NHAA*6,-200000);
        natom.resize(NHAA,0);
        soften.resize(NHAA,0);
        hlig.resize(NHAA,-1);
        ih.resize(N,-1);
        seq.resize(N,NONE);
        chainof.resize(NHAA,-1);
        Mehp.resize(NHAA*NHAA,0);
        Vehp.resize(NHAA*NHAA,0);
        changed.resize(NHAA*NHAA,0);
        Ligand *lg;

        for (int i1=0,k1=0,l1=0;i1<NC();++i1) {

            for (int j1=0;j1<p->Chain(i1)->numLigands();++j1,++k1) {
                lg=p->Chain(i1)->memberLigand(j1);

                if (htyp(seq[k1]=lg->OLC()) ==-1) continue;
                
                hlig[l1] = lg->UniqueId(); // for each hydrophobic store ligand index
                 
                 
                 prf::cout << "HydrophobicityPy::init "<< Groups::mapOLC2Char(lg->OLC()) << "," << lg->UniqueId()+1 << " - "<<l1<<"\n";

                
                for (int i2=0,k2=0,l2=0;i2<NC();++i2) {
                    for (int j2=0;j2<p->Chain(i2)->numLigands();++j2,++k2) {
                        if (htyp(seq[k2]=p->Chain(i2)->memberLigand(j2)->OLC())==-1) continue;

                        tmpx=hpstr[htyp(seq[k1])] + hpstr[htyp(p->ligand(k2)->OLC())];

                        if (i1==i2&&(k2==(k1+1)||k2==(k1-1))) tmpx*=strnn;

                        if (i1==i2&&(k2==(k1+2)||k2==(k1-2))) tmpx*=strnnn;

                        strength[l1*NHAA+l2++]=tmpx;
                    }
                }

                chainof[l1]=i1;

                if (seq[k1]==P || seq[k1]==F || seq[k1]==Y
                    || seq[k1]==W || seq[k1]==DPR) soften[l1]=1;
                else if (seq[k1]==L || seq[k1]==I || seq[k1]==V) soften[l1]=2;

                ih[k1]=l1++;
            }
        }

        blog(50)<<"Hydrophobic Ligands list ...\n";

        for (int i=0;i<p->NumberOfLigands();++i) {
            if (ih[i]>=0) {
                blog<<ih[i]<<"  "<<Groups::mapOLC2Char(seq[i])<<" on chain "
                <<chainof[ih[i]]<<" index "<<p->ligand(i)->SeqSerial()
                <<"\n";
            }
        }

        blog<<"Pairwise hydrophobicity strength matrix ...\n";

        for (int i=0;i<NHAA;++i) {
            for (int j=0;j<NHAA;++j) {
                blog<<strength[i*NHAA+j]<<"  ";
            }

            blog<<"\n";
        }

        for (int ilg=0,n=0;ilg<p->NumberOfLigands();++ilg) {
            if ((n=ih[ilg]) <0) continue;

            for (int i=0;i<6;++i) {
                Atom *at1=p->ligand(ilg)->labeled_atom(
                              hatms[htyp(seq[ilg])][i]);

                if (at1!=NULL) {atom[6*n+i]=at1->UniqueId();natom[n]++;}
            }
        }

        double max_at_dist=(14+sqrt(a2));

        max_at_dist2=max_at_dist*max_at_dist;
        max_at_dist=AtomCoordinates::boxL()-max_at_dist;
        max_cmp_dst2=max_at_dist*max_at_dist;
        initialized=true;
    }

    void Hydrophobicity::DisplayMatrix()
    {
        for (int i=0;i<NHAA;++i) {
            for (int j=0;j<NHAA;++j) prf::cout <<Mehp[NHAA*i+j]<<"\t";

            prf::cout <<"\n";
        }

        for (int i=0;i<NHAA;++i) {
            for (int j=0;j<NHAA;++j) prf::cout <<Vehp[NHAA*i+j]<<"\t";

            prf::cout <<"\n";
        }
    }

    double Hydrophobicity::evaluate()
    {
        double eterm=0;
        vval=delv=0;
        nchanges=0;

        for (int i=0;i<NHAA;++i) {
            for (int j=i+1;j<NHAA;++j) {
                eterm=hp_pair(i,j,chainof[i]==chainof[j]?
                              AtomCoordinates::s2:
                              AtomCoordinates::dist2);
                vval+=eterm;
                Mehp[NHAA*i+j]=Mehp[NHAA*j+i]=eterm;
            }
        }

        return vval*=sclfct;
    }

    double Hydrophobicity::gradientXYZ(std::valarray<double> &gx)
    {
        double eterm=0;
        vval=delv=0;
        nchanges=0;
        gx=0;
        for (int i=0;i<NHAA;++i) {
            for (int j=i+1;j<NHAA;++j) {
                eterm=dhp_pair(i,j,chainof[i]==chainof[j],gx);
                vval+=eterm;
                Mehp[NHAA*i+j]=Mehp[NHAA*j+i]=eterm;
            }
        }

        return vval*=sclfct;

    }

    void Hydrophobicity::rangeEstimate(double &x1, double &x2)
    {

        x1=-1.2*sclfct*max(NHAA,1);

        x2=0;
        //mere guess.
    }

    double Hydrophobicity::deltaE(Update *updt)
    {
        delv=0;
        nchanges=0;

        if (updt->sidechain_update()) {
            int ind=ih[updt->change(0).info.group];

            if (ind<0 || atom[6*ind+natom[ind]-1]<updt->begin_atom() ||
                atom[6*ind]>=updt->end_atom()) return 0;
        }

        double eterm=0,de=0;

        updt->n_residue_rigid_ranges();
        vector<pair<int,int> > *rng=updt->residue_rigid_ranges();

        for (int it=0; it<updt->n_residue_rigid_ranges();++it) {
            int ires1=(*rng)[it].first, ires2=(*rng)[it].second+1;

            for (int jres=0,prng=0;jres<N;) {
                int ihres=0,jhres=0;

                if ((jhres=ih[jres])<0) {++jres;continue;} // residue j is not hydrophobic, skip.

                if (prng<=it && jres>=(*rng)[prng].first) {
                    jres=(*rng)[prng++].second+1;
                    continue;
                }

                for (int ires=ires1;ires<ires2;++ires) {
                    if ((ihres=ih[ires])<0) continue;

                    if (chainof[ihres]==chainof[jhres])
                        eterm=hp_pair(ihres,jhres,AtomCoordinates::s2);
                    else eterm=hp_pair(ihres,jhres,AtomCoordinates::dist2);

                    delv+=(de=eterm-Mehp[ihres*NHAA+jhres]);

                    if (de!=0) {
                        Vehp[changed[nchanges++]=ihres*NHAA+jhres]=eterm;
                        Vehp[changed[nchanges++]=jhres*NHAA+ihres]=eterm;
                    }
                }

                ++jres;
            }
        }

        return (delv*=sclfct);
    }

    void Hydrophobicity::Accept(Update *updt)
    {
        for (;nchanges;--nchanges)
            Mehp[changed[nchanges-1]]=Vehp[changed[nchanges-1]];

        vval+=delv;
    }

    double Hydrophobicity::hp_pair(int i, int j, double(*distf)(int,int))
    {
        /* i and j are indices for NHAA sized arrays. */
        double ans=-strength[i*NHAA+j]*contact_frac(i,j,distf);

        if (ans>0) prf::cerr <<"hp_pair: "<<i<<"\t"<<j<<" : "<<ans<<"\n";
        
        
        if (ans != 0.0 && hlig[i] != -1 &&  hlig[j] != -1){
        	/*char olc1 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[i]))->OLC() );
	        int  id1  = ((AminoAcid*) p->amino_acid(hlig[i]))->UniqueId();
	        char olc2 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[j]))->OLC() );
	        int  id2  = ((AminoAcid*) p->amino_acid(hlig[j]))->UniqueId();*/
        	
	        if ( ( hlig[i]==focalRes && hlig[j]==focalRes2)
	        		|| ( hlig[i]==focalRes2 && hlig[j]==focalRes) ){
	        	ans *= focalResScale;
	        }
	        /*
        	prf::cout << "HydrophobicityPi::hp_pair Res:"<<olc1<<id1+1 <<",Res:"<<olc2<< id2+1 << ": " << ans << "";
        	if ( hlig[i]==focalRes || hlig[j]==focalRes){
        		prf::cout << " ***\n";
        		
	        	//if ( hlig[j] != -1 && ((AminoAcid*) p->amino_acid(hlig[j]))->OLC() == F && hlig[j]==51){
	        	//	ans *= 1.0;
	        	//prf::cout << "Pi:" << hlig[i] << "," << hlig[j] << "\n";
	        
	        	//}
        	}else{
	        	prf::cout << "\n";
	        }
	        */
        }
        
        
        return ans;
    }

    double Hydrophobicity::hp_contact_frac(int iaa, int jaa,
                                           double(*distf)(int,int))
    {
        /* iaa and jaa are amino acid indices relative to the
        entire system. This function is intended only for use
        in connection with hydrophobicity related observables.*/
        int i=ih[iaa],j=ih[jaa];

        if (i<0||j<0) return 0.0;

        return contact_frac(i,j,distf);
    }

    double Hydrophobicity::contact_frac(int i, int j,
                                        double(*distf)(int,int))
    {
        /* i and j are indices for NHAA sized arrays. */
        int ni,nj,k,l,m,n;
        double tmp,sum=0;
        ni=natom[i];
        nj=natom[j];

        for (n=0;n<ni+nj;n++) r2min[n]=1e10;

        for (n=0;n<ni;n++) {
            k=atom[i*6+n];

            for (m=0;m<nj;m++) {
                l=atom[j*6+m];
                tmp=distf(k,l); // distance between atoms k,l

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

        double gammaij=1;

        if ((soften[i] == soften[j]) and(soften[i]!=0)) gammaij=0.75;

        return std::min(1.0,sum/(gammaij*(ni+nj)));
    }
}

double Hydrophobicity::dhp_pair(int i, int j, bool perdhint,
                                std::valarray<double> &gx)
{
    double frac=0,pstr=0;
    int ni,nj,k,l,m,n;
    double tmp,sum=0;
    ni=natom[i];
    nj=natom[j];
    pstr=-strength[i*NHAA+j]*sclfct;
    for (n=0;n<ni+nj;++n) r2min[n]=1e10;

    for (n=0;n<ni;++n) {
        k=atom[i*6+n];

        for (m=0;m<nj;++m) {
            l=atom[j*6+m];
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
        else sum+= (b2-r2min[n]) /d2;
    }

    double gammaij=1;

    if ((soften[i] == soften[j]) and(soften[i]!=0)) gammaij=0.75;

    frac=(gammaij*(ni+nj));
    if (sum<=0) return 0;
    else if (sum<frac) {
        pstr/=(frac);
        for (n=0;n<ni+nj;++n) {
            if (r2min[n]>b2 or r2min[n]<a2) continue;
            if (n<ni) {
                m=atom[i*6+n];
                l=atom[6*j+r2minpartner[n]];
            } else {
                m=atom[6*j+n-ni];
                l=atom[6*i+r2minpartner[n]];
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
