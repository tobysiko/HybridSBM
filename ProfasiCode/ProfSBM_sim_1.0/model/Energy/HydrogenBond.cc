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

#include "HydrogenBond.hh"
using std::vector;
using std::valarray;
using std::pair;
using std::max;

namespace prf
{
    ///////////////HydrogenBond Base Class ////////////////////////
    HydrogenBond::HydrogenBond()
    {
        sighb=2.0;
        cuthb=4.5;
        cut2=cuthb*cuthb;
        POWA=0.5;
        POWB=0.5;
        epshb1=3.0;
        epshb2=2.3;
        bhb=-30*(pow(sighb/cuthb,10.0)-pow(sighb/cuthb,12.0))/cut2;
        ahb=-(5*pow(sighb/cuthb,12.0)-6*pow(sighb/cuthb,10.0))-bhb*cut2;
        sighb2=sighb*sighb;
        halfpowa=0.5*POWA;
        halfpowb=0.5*POWB;
        ntdonsupr=1.0;
        ntaccsupr=1.0;
        ctdonsupr=1.0;
        ctaccsupr=1.0;
        cdon=pow(1.0/AminoAcid::bNH,POWA);
        cacc=pow(1.0/AminoAcid::bCO,POWB);
        csacc=pow(1.0/1.25,POWB);
    }

    HydrogenBond::~HydrogenBond(){}

    void HydrogenBond::set_population(Population *p)
    {
        ndon=nacc=0;
        double termdonsupr=0.5, termaccsupr=0.5; //Ugly and temporary

        for (int ich=0;ich<p->NumberOfChains();++ich) {
            for (int i=0;i<p->Chain(ich)->numAminoAcids();++i) {
                AminoAcid * ac=p->Chain(ich)->AA(i);

                if (!ac->hasNTerminal()) ndon+=ac->n_donors();

                if (!ac->hasCTerminal()) nacc+=ac->n_acceptors();
            }

            if (p->Chain(ich)->NtermLigand()!=NULL) {
                ndon+=p->Chain(ich)->NtermLigand()->n_donors();
                nacc+=p->Chain(ich)->NtermLigand()->n_acceptors();
            }

            if (p->Chain(ich)->CtermLigand()!=NULL) {
                ndon+=p->Chain(ich)->CtermLigand()->n_donors();
                nacc+=p->Chain(ich)->CtermLigand()->n_acceptors();
            }
        }

        donor.resize(ndon);

        acceptor.resize(nacc);

        donstart.resize(p->NumberOfLigands(),-1);
        accstart.resize(p->NumberOfLigands(),-1);
        donend.resize(p->NumberOfLigands(),-1);
        accend.resize(p->NumberOfLigands(),-1);
        int idon=0,iacc=0;
        Dipole d;

        for (int ich=0;ich<p->NumberOfChains();++ich) {
            for (int ilg=0;ilg<p->Chain(ich)->numLigands();++ilg) {
                Ligand * lg=p->Chain(ich)->memberLigand(ilg);
                donstart[lg->UniqueId()]=idon;
                accstart[lg->UniqueId()]=iacc;

                if (lg->isAA()) {
                    int i=lg->UniqueId()-p->Chain(ich)->AA(0)->UniqueId();
                    AminoAcid * acd=p->Chain(ich)->AA(i);

                    if (!acd->hasNTerminal()) {
                        for (int j=0;j<acd->n_donors();++j) {
                            acd->Donor(j,d);
                            d.ligand(acd->UniqueId());
                            d.chain(acd->LocatedOn());
                            d.SetStrength(cdon);
                            donor[idon++]=d;
                        }
                    }

                    if (!acd->hasCTerminal()) {
                        for (int j=0;j<acd->n_acceptors();++j) {
                            acd->Acceptor(j,d);
                            d.ligand(acd->UniqueId());
                            d.chain(acd->LocatedOn());
                            d.SetStrength(cacc);
                            acceptor[iacc++]=d;
                        }
                    }
                } else {
                    for (int id=0;id<lg->n_donors();++id) {
                        lg->Donor(id,d);
                        d.ligand(lg->UniqueId());
                        d.chain(lg->LocatedOn());
                        d.SetStrength(termdonsupr*cdon);
                        donor[idon++]=d;
                    }

                    for (int id=0;id<lg->n_acceptors();++id) {
                        lg->Acceptor(id,d);
                        d.ligand(lg->UniqueId());
                        d.chain(lg->LocatedOn());
                        d.SetStrength(termaccsupr*cacc);
                        acceptor[iacc++]=d;
                    }
                }

                donend[lg->UniqueId()]=idon;

                accend[lg->UniqueId()]=iacc;
            }

            if (p->Chain(ich)->NtermLigand()!=NULL) {
                int ist=p->Chain(ich)->AA(0)->UniqueId(),ind=0;
                ind=donend[ist];
                ist=donstart[ist];

                for (int id=ist;id<ind;++id) donor[id].Reduce(termdonsupr);
            }

            if (p->Chain(ich)->CtermLigand()!=NULL) {
                int ist=p->Chain(ich)->AA(p->Chain(ich)->numAminoAcids()-1)->UniqueId(),ind=0;
                ind=accend[ist];
                ist=accstart[ist];

                for (int ia=ist;ia<ind;++ia) acceptor[ia].Reduce(termaccsupr);
            }

        }
    }

    double HydrogenBond::dhb_pair(Dipole &don, Dipole &acc, int excl1,int excl2,
                                  std::valarray<double> &gx)
    {
        double r,r2,r4,r6;
        double ca,cb,scacb,Vr,dVr;
        Vector3 dv,d12,a12,dirn;
        int d1=don.atom1(),d2=don.atom2(),a1=acc.atom1(),a2=acc.atom2();

        if (don.chain()!=acc.chain()) {
            dv=AtomCoordinates::sep(a1,d2);
            d12=AtomCoordinates::sep(d2,d1);
            a12=AtomCoordinates::sep(a2,a1);
        } else {
            if ((acc.ligand()-don.ligand())>=excl1 &&
                (acc.ligand()-don.ligand())<=excl2) return 0;

            dv=AtomCoordinates::diff(a1,d2);
            d12=AtomCoordinates::diff(d2,d1);
            a12=AtomCoordinates::diff(a2,a1);
        }

        if ((r2=dv.mag2()) > cut2) return 0;
        if ((ca=d12.dot(dv))<=0) return 0;
        if ((cb=a12.dot(dv))<=0) return 0;

        scacb=sqrt(ca*cb);
        r=sqrt(r2);
        r6=sighb2/r2;
        r4=r6*r6;
        r6=r6*r4;
        Vr=(r6*(5*r6-6*r4)+ahb+bhb*r2)/r;
        dVr=(2*bhb*r-Vr-(60/r)*(r6*(r6-r4)))/r2;
        Vr*=(don.strength()*acc.strength());
        dVr*=(don.strength()*acc.strength());

        dirn=dVr*scacb*dv;
        gx[3*a1]+=dirn.x();
        gx[3*a1+1]+=dirn.y();
        gx[3*a1+2]+=dirn.z();
        gx[3*d2]-=dirn.x();
        gx[3*d2+1]-=dirn.y();
        gx[3*d2+2]-=dirn.z();

        dVr=Vr/2.0/scacb;

        dirn=cb*dVr*d12+ca*dVr*(a12-dv);
        gx[3*a1]+=dirn.x();
        gx[3*a1+1]+=dirn.y();
        gx[3*a1+2]+=dirn.z();

        dirn=ca*dVr*dv;
        gx[3*a2]+=dirn.x();
        gx[3*a2+1]+=dirn.y();
        gx[3*a2+2]+=dirn.z();

        dirn=-cb*dVr*dv;
        gx[3*d1]+=dirn.x();
        gx[3*d1+1]+=dirn.y();
        gx[3*d1+2]+=dirn.z();

        dirn=cb*dVr*(dv-d12)-ca*dVr*a12;
        gx[3*a1]+=dirn.x();
        gx[3*a1+1]+=dirn.y();
        gx[3*a1+2]+=dirn.z();

        return Vr*scacb;
    }

    double HydrogenBond::hb_pair(Dipole &don, Dipole &acc, int excl1,int excl2)
    {
        double r2,r4,r6;
        double ca,cb;
        Vector3 dv;
        int d1=don.atom1(),d2=don.atom2(),a1=acc.atom1(),a2=acc.atom2();

        if (don.chain()!=acc.chain()) {
            dv=AtomCoordinates::sep(a1,d2);

            if ((r2=dv.mag2()) > cut2) return 0;

            ca=AtomCoordinates::sep(d2,d1).dot(dv);

            cb=AtomCoordinates::sep(a2,a1).dot(dv);
        } else {
            if ((acc.ligand()-don.ligand())>=excl1 &&
                (acc.ligand()-don.ligand())<=excl2) return 0;

            dv=AtomCoordinates::diff(a1,d2);

            if ((r2=dv.mag2()) > cut2) return 0;

            ca=AtomCoordinates::diff(d2,d1).dot(dv);

            cb=AtomCoordinates::diff(a2,a1).dot(dv);
        }

        if (POWA>0 &&ca<0) return 0;

        if (POWB>0 &&cb<0) return 0;

        r6=sighb2/r2;

        r4=r6*r6;

        r6=r6*r4;

        return don.strength()*acc.strength()*
               pow(ca*ca/r2,halfpowa)*pow(cb*cb/r2,halfpowb)*
               (r6*(5*r6-6*r4)+ahb+bhb*r2);
    }

    ///////////////HBMM////////////////////////////////
    HBMM::HBMM()
    {
        Name("HBMM");
        grdtyp=2;
    }

    HBMM::~HBMM() {}

    void HBMM::init()
    {
        if (initialized) return;
        HydrogenBond::set_population(p);
        Mmm.resize(ndon*nacc,0);
        Vmm.resize(ndon*nacc,0);
        changed.resize(ndon*nacc,0);
        initialized=true;
    }

    void HBMM::rangeEstimate(double &x1, double &x2)
    {
        int npro=0,nendgrp=0,Naa=0;

        for (int j=0;j<NC();++j) {
            for (int i=0;i<p->Chain(j)->numAminoAcids();++i) {
                if (p->Chain(j)->memberAA(i)->OLC()==P) ++npro;

                if (p->Chain(j)->memberAA(i)->OLC()==DPR) ++npro;

                ++Naa;
            }

            if (p->Chain(j)->NtermLigand()!=NULL) ++nendgrp;

            if (p->Chain(j)->CtermLigand()!=NULL) ++nendgrp;
        }

        x1=-2.3*(Naa-0.5*(npro+2*NC()-nendgrp));

        x2=-0.01*x1;
        //guess.
    }

    double HBMM::evaluate()
    {
        double eterm=0;
        vval=delv=0;

        for (int i=0;i<ndon;++i) {
            for (int j=0;j<nacc;++j) {
                eterm=hb_pair(donor[i],acceptor[j],-2,0);
                vval+=(Mmm[nacc*i+j]=eterm);
            }
        }

        return vval*=epshb1;
    }

    double HBMM::gradientXYZ(std::valarray<double> &gx)
    {
        double eterm=0;
        vval=delv=0;
        gx=0;

        for (int i=0;i<ndon;++i) {
            for (int j=0;j<nacc;++j) {
                eterm=dhb_pair(donor[i],acceptor[j],-2,0,gx);
                vval+=eterm;
            }
        }

        gx*=epshb1;

        return vval*=epshb1;
    }

    // iaa and jaa are global ligand indices
    double HBMM::InterLg(int iaa, int jaa)
    {
        double etot=0;

        for (int i=donstart[iaa];i<donend[iaa];++i) {
            for (int j=accstart[jaa];j<accend[jaa];++j) {
                etot+=hb_pair(donor[i],acceptor[j],-2,0);
            }
        }

        return etot*=epshb1;
    }

    double HBMM::InterAA(int iaa, int jaa)
    {
        return InterLg(p->amino_acid(iaa)->UniqueId(),
                       p->amino_acid(jaa)->UniqueId());
    }

    double HBMM::InterAA(int ich, int jch, int iaa, int jaa)
    {
        return InterLg(p->Chain(ich)->AA(iaa)->UniqueId(),
                       p->Chain(jch)->AA(jaa)->UniqueId());
    }

    double HBMM::InterLg(int ich, int jch, int iaa, int jaa)
    {
        return InterLg(p->Chain(ich)->memberLigand(iaa)->UniqueId(),
                       p->Chain(jch)->memberLigand(jaa)->UniqueId());
    }

    double HBMM::InterChain(int ich, int jch)
    {
        int l1=p->Chain(ich)->first_ligand()->UniqueId();
        int l2=p->Chain(ich)->last_ligand()->UniqueId();
        int l3=p->Chain(jch)->first_ligand()->UniqueId();
        int l4=p->Chain(jch)->last_ligand()->UniqueId();
        double etot=0;

        for (int i=donstart[l1];i<donend[l2];++i) {
            for (int j=accstart[l3];j<accend[l4];++j) {
                etot+=hb_pair(donor[i],acceptor[j],-2,0);
            }
        }

        return etot*=epshb1;
    }

    double HBMM::interhb()
    {
        double e=0;

        for (int i=0;i<NC();i++)
            for (int j=0;j<NC();j++) {
                if (i!=j) e+=InterChain(i,j);
            }

        return e;
    }

    void HBMM::Accept(Update *updt)
    {
        vval+=delv;

        for (; nchanges; --nchanges)
            Mmm[changed[nchanges-1]]=Vmm[changed[nchanges-1]];
    }

    double HBMM::deltaE(Update *updt)
    {
        delv=0;
        nchanges=0;

        if (updt->sidechain_update() &&
            p->ligand((updt->change(0).info.group))->isAA())
            return delv;

        double eterm=0,de=0;

        int numsets=updt->n_residue_rigid_ranges();

        vector<pair<int,int> > *rng=updt->residue_rigid_ranges();

        for (int i=0;i<numsets;++i) {
            int i1=(*rng)[i].first,i2=(*rng)[i].second;

            for (int jacc=0,prng=0;jacc<nacc;) {
                if (prng<=i && jacc>=accstart[(*rng)[prng].first]) {
                    jacc=accend[(*rng)[prng++].second];
                    continue;
                }

                for (int idon=donstart[i1];idon<donend[i2];++idon) {
                    eterm=hb_pair(donor[idon],acceptor[jacc],-2,0);
                    delv+=(de=eterm-Mmm[nacc*idon+jacc]);

                    if (de!=0) {
                        Vmm[nacc*idon+jacc]=eterm;
                        changed[nchanges++]=nacc*idon+jacc;
                    }
                }

                ++jacc;
            }

            for (int jdon=0,prng=0;jdon<ndon;) {
                if (prng<=i && jdon>=donstart[(*rng)[prng].first]) {
                    jdon=donend[(*rng)[prng++].second];
                    continue;
                }

                for (int iacc=accstart[i1];iacc<accend[i2];++iacc) {
                    eterm=hb_pair(donor[jdon],acceptor[iacc],-2,0);
                    delv+=(de=eterm-Mmm[nacc*jdon+iacc]);

                    if (de!=0) {
                        Vmm[nacc*jdon+iacc]=eterm;
                        changed[nchanges++]=nacc*jdon+iacc;
                    }
                }

                ++jdon;
            }
        }

        return delv*=epshb1;
    }


    ////////////////////HBMS//////////////////////////
    HBMS::HBMS()
    {
        Name("HBMS");
        grdtyp=2;
    }

    HBMS::~HBMS() {}

    void HBMS::rangeEstimate(double &x1, double &x2)
    {
        int ncontribs=0;

        for (int i=0;i<p->NumberOfLigands();++i) {
            if (participates(p->ligand(i)->OLC())) ++ncontribs;
        }

        for (int i=0;i<p->NumberOfChains();++i) {
            if (p->Chain(i)->AA(0)->hasNTerminal()) ++ncontribs;
            if (p->Chain(i)->AA(p->Chain(i)->numAminoAcids()-1)->hasCTerminal())
                ++ncontribs;
        }

        x1=-1.5*max(ncontribs,1);
        x2=-0.01*x1;
        /*
          A perfectly satisfied main chain side chain HB will contribute
          about 2 units. The 1.5 says that we don't expect all of the
          charged side chains to be involved in HB with the backbone
          simultaneously. Each charged residue contains more than one
          dipoles, but again, it is not possible (or very unlikely)
          that they all participate in HB with the backbone simultaneously.
          So, we count the number of possible contributors: charged
          side chains and chain ends. Then we count about 2 units for each,
          and throw in a 0.75 to account for the fact that some side chains
          will normally miss an HBMS. The estimate is linear in the number
          of contributors.
          */
    }

    bool HBMS::participates(prf::OneLetterCode cd)
    {
        return cd==D or cd==E or cd==SUC or cd==K or cd==R;
    }

    void HBMS::init()
    {
        if (initialized) return;
        HydrogenBond::set_population(p);
        nscdon=nscacc=0;

        for (int i=0;i<p->NumberOfLigands();++i) {
            Ligand *lg=p->ligand(i);

            switch (lg->OLC()) {
                case D:
                case E:
                case SUC:
                    nscacc+=2;
                    break;
                case K:
                    nscdon+=3;
                    break;
                case R:
                    nscdon+=5;
                    break;
                default:
                    break;
            };
        }

        for (int i=0;i<p->NumberOfChains();++i) {
            if (p->Chain(i)->AA(0)->hasNTerminal()) nscdon+=3;

            if (p->Chain(i)->AA(p->Chain(i)->numAminoAcids()-1)->hasCTerminal()) nscacc+=2;
        }

        scdonor.resize(nscdon);

        scacceptor.resize(nscacc);
        scdonbeg.resize(p->NumberOfLigands(),0);
        scdonend.resize(p->NumberOfLigands(),0);
        scaccbeg.resize(p->NumberOfLigands(),0);
        scaccend.resize(p->NumberOfLigands(),0);
        int idon=0,iacc=0;

        for (int ich=0;ich<p->NumberOfChains();++ich) {
            for (int ilg=0;ilg<p->Chain(ich)->numLigands();++ilg) {
                Ligand * ac = p->Chain(ich)->memberLigand(ilg);
                int i=ac->UniqueId();
                scdonbeg[i]=idon;
                scaccbeg[i]=iacc;
                Dipole d;

                if (ac->isAA() && ilg==0) {
                    d=Dipole(ac->at(" N  ").UniqueId(),
                             ac->at("1H  ").UniqueId(),cdon);
                    d.ligand(ac->UniqueId());
                    d.chain(ac->LocatedOn());
                    scdonor[idon++]=d;
                    d=Dipole(ac->at(" N  ").UniqueId(),
                             ac->at("2H  ").UniqueId(),cdon);
                    d.ligand(ac->UniqueId());
                    d.chain(ac->LocatedOn());
                    scdonor[idon++]=d;
                    d=Dipole(ac->at(" N  ").UniqueId(),
                             ac->at("3H  ").UniqueId(),cdon);
                    d.ligand(ac->UniqueId());
                    d.chain(ac->LocatedOn());
                    scdonor[idon++]=d;
                }

                if (ac->isAA() && ilg==(p->Chain(ich)->numLigands()-1)) {
                    d=Dipole(ac->at(" O  ").UniqueId(),
                             ac->at(" C  ").UniqueId(),cacc);
                    d.ligand(ac->UniqueId());
                    d.chain(ac->LocatedOn());
                    scacceptor[iacc++]=d;
                    d=Dipole(ac->at(" OXT").UniqueId(),
                             ac->at(" C  ").UniqueId(),cacc);
                    d.ligand(ac->UniqueId());
                    d.chain(ac->LocatedOn());
                    scacceptor[iacc++]=d;
                }

                switch (ac->OLC()) {
                    case D: {
                        d=Dipole(ac->at(" OD1").UniqueId(),
                                 ac->at(" CG ").UniqueId(),csacc);
                        d.ligand(p->ligand(i)->UniqueId());
                        d.chain(p->ligand(i)->LocatedOn());
                        scacceptor[iacc++]=d;
                        d=Dipole(ac->at(" OD2").UniqueId(),
                                 ac->at(" CG ").UniqueId(),csacc);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scacceptor[iacc++]=d;
                        break;
                    }

                    case E: {
                        d=Dipole(ac->at(" OE1").UniqueId(),
                                 ac->at(" CD ").UniqueId(),csacc);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scacceptor[iacc++]=d;
                        d=Dipole(ac->at(" OE2").UniqueId(),
                                 ac->at(" CD ").UniqueId(),csacc);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scacceptor[iacc++]=d;
                        break;
                    }

                    case SUC: {
                        d=Dipole(ac->at(" O11").UniqueId(),
                                 ac->at(" CO2").UniqueId(),csacc);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scacceptor[iacc++]=d;
                        d=Dipole(ac->at(" O12").UniqueId(),
                                 ac->at(" CO2").UniqueId(),csacc);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scacceptor[iacc++]=d;
                        break;
                    }

                    case K: {
                        d=Dipole(ac->at(" NZ ").UniqueId(),
                                 ac->at("1HZ ").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        d=Dipole(ac->at(" NZ ").UniqueId(),
                                 ac->at("2HZ ").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        d=Dipole(ac->at(" NZ ").UniqueId(),
                                 ac->at("3HZ ").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        break;

                    }

                    case R: {
                        d=Dipole(ac->at(" NE ").UniqueId(),
                                 ac->at(" HE ").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        d=Dipole(ac->at(" NH1").UniqueId(),
                                 ac->at("1HH1").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        d=Dipole(ac->at(" NH1").UniqueId(),
                                 ac->at("2HH1").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        d=Dipole(ac->at(" NH2").UniqueId(),
                                 ac->at("1HH2").UniqueId(),cdon);
                        d.ligand(ac->UniqueId());
                        d.chain(ac->LocatedOn());
                        scdonor[idon++]=d;
                        d=Dipole(ac->at(" NH2").UniqueId(),
                                 ac->at("2HH2").UniqueId(),cdon);
                        d.chain(ac->LocatedOn());
                        d.ligand(ac->UniqueId());
                        scdonor[idon++]=d;
                        break;
                    }

                    default:
                        break;
                };

                scdonend[i]=idon;

                scaccend[i]=iacc;
            }
        }

        Vms.resize(nscacc*ndon+nscdon*nacc,0);

        Mms.resize(nscacc*ndon+nscdon*nacc,0);
        changed.resize(nscacc*ndon+nscdon*nacc,0);
        initialized=true;
    }

    double HBMS::gradientXYZ(std::valarray<double> &gx)
    {
        vval=delv=0;
        double eterm=0;

        for (int i=0;i<nscacc;++i) {
            for (int j=0;j<ndon;++j) {
                eterm=dhb_pair(donor[j],scacceptor[i],1,-1,gx);
                vval+=eterm;
            }
        }

        for (int i=0;i<nscdon;++i) {
            for (int j=0;j<nacc;++j) {
                eterm=dhb_pair(scdonor[i],acceptor[j],1,-1,gx);
                vval+=eterm;
            }
        }

        gx*=epshb2;
        return vval*=epshb2;
    }

    double HBMS::evaluate()
    {
        vval=delv=0;
        nchanges=0;
        double eterm=0;

        for (int i=0;i<nscacc;++i) {
            for (int j=0;j<ndon;++j) {
                eterm=hb_pair(donor[j],scacceptor[i],1,-1);
                vval+=(Mms[ndon*i+j]=eterm);
            }
        }

        for (int i=0;i<nscdon;++i) {
            for (int j=0;j<nacc;++j) {
                eterm=hb_pair(scdonor[i],acceptor[j],1,-1);
                vval+=(Mms[ndon*nscacc+nacc*i+j]=eterm);
            }
        }

        return vval*=epshb2;
    }

    void HBMS::Accept(Update *updt)
    {
        vval+=delv;

        for (; nchanges; --nchanges)
            Mms[changed[nchanges-1]]=Vms[changed[nchanges-1]];
    }

    double HBMS::deltaE(Update *updt)
    {
        Ligand *lg=p->ligand(updt->change(0).info.group);
        OneLetterCode o=lg->OLC();
        delv=0;
        nchanges=0;

        if (updt->sidechain_update() && (not participates(o))) return delv;

        double eterm=0,de=0;

        int idum=0,numsets=updt->n_residue_rigid_ranges();

        vector<pair<int,int> > *rng=updt->residue_rigid_ranges();

        for (int i=0;i<numsets;++i) {
            int i1=(*rng)[i].first,i2=(*rng)[i].second;

            for (int id=0,prng=0;id<ndon;) {
                if (prng<i&& id>=donstart[(*rng)[prng].first]) {
                    id=donend[(*rng)[prng++].second];
                    continue;
                }

                for (int ja=scaccbeg[i1];ja<scaccend[i2];++ja) {
                    eterm=hb_pair(donor[id],scacceptor[ja],1,-1);
                    delv+=(de=eterm-Mms[idum=ndon*ja+id]);

                    if (de!=0) {
                        Vms[idum]=eterm;
                        changed[nchanges++]=idum;
                    }
                }

                ++id;
            }

            for (int ja=0,prng=0;ja<nacc;) {
                if (prng<i && ja>=accstart[(*rng)[prng].first]) {
                    ja=accend[(*rng)[prng++].second];
                    continue;
                }

                for (int id=scdonbeg[i1];id<scdonend[i2];++id) {
                    eterm=hb_pair(scdonor[id],acceptor[ja],1,-1);
                    delv+=(de=eterm-Mms[idum=nscacc*ndon+nacc*id+ja]);

                    if (de!=0) {
                        Vms[idum]=eterm;
                        changed[nchanges++]=idum;
                    }
                }

                ++ja;
            }

            if (updt->sidechain_update() && lg->isAA()) continue;

            for (int ja=0,prng=0;ja<nscacc;) {
                if (prng<=i && ja>=scaccbeg[(*rng)[prng].first]) {
                    ja=scaccend[(*rng)[prng++].second];
                    continue;
                }

                for (int id=donstart[i1];id<donend[i2];++id) {
                    eterm=hb_pair(donor[id],scacceptor[ja],1,-1);
                    delv+=(de=eterm-Mms[idum=ndon*ja+id]);

                    if (de!=0) {
                        Vms[idum]=eterm;
                        changed[nchanges++]=idum;
                    }
                }

                ++ja;
            }

            for (int id=0,prng=0;id<nscdon;) {
                if (prng<=i && id>=scdonbeg[(*rng)[prng].first]) {
                    id=scdonend[(*rng)[prng++].second];
                    continue;
                }

                for (int ja=accstart[i1];ja<accend[i2];++ja) {
                    eterm=hb_pair(scdonor[id],acceptor[ja],1,-1);
                    delv+=(de=eterm-Mms[idum=ndon*nscacc+nacc*id+ja]);

                    if (de!=0) {
                        Vms[idum]=eterm;
                        changed[nchanges++]=idum;
                    }
                }

                ++id;
            }
        }

        return delv*=epshb2;
    }
}
