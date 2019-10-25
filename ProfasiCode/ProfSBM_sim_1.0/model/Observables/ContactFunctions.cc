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

#include "ContactFunctions.hh"

namespace prf
{
    ContactFunction::ContactFunction() : minseqsep(2), issym(true) {}

    ContactFunction::~ContactFunction() {}

    int ContactFunction::init(Population *p) {return 0;}

    void ContactFunction::set_cutoff(double xx) {}

    double ContactFunction::get_cutoff() const {return 0;}

    bool ContactFunction::operator()(int arg1,int arg2) const {return true;}
    
    double ContactFunction::operator()(int arg1,int arg2, bool dummy) const {return -2.0;} // NEW!! -TS
    
    bool ContactFunction::operator()(const Contact & c) const
    {
        return operator()(c.first,c.second);
    }

    CaContact::CaContact()
    {
        setCut(6.0);
    }

    CaContact::CaContact(double dtc)
    {
        setCut(dtc);
    }

    void CaContact::set_cutoff(double xx)
    {
        setCut(xx);
    }

    double CaContact::get_cutoff() const {return distcut;}

    int CaContact::init(Population *p)
    {
        indca.resize(p->NumberOfLigands(),-1);

        for (size_t i=0;i<indca.size();++i) {
            Atom * atm;

            if ((atm=p->ligand(i)->labeled_atom(" CA "))!=NULL)
                indca[i]=atm->UniqueId();
        }

        distcut2=distcut*distcut;

        return 1;
    }

    bool CaContact::operator()(const Contact &c) const
    {
        return operator()(c.first,c.second);
    }

    bool CaContact::operator()(int i1,int i2) const
    {
        if (abs(i1-i2) <minseqsep) return false;

        if (indca[i1]==-1 || indca[i2]==-1) return false;

        return (AtomCoordinates::dist2(indca[i1],indca[i2]) <distcut2);
    }
    
    double CaContact::operator()(int i1,int i2, bool dummy) const //  NEW! -TS
        {
            if (abs(i1-i2) <minseqsep) return -1.0;

            if (indca[i1]==-1 || indca[i2]==-1) return -1.0;
            
            double distance = AtomCoordinates::dist2(indca[i1],indca[i2]);
            //std::cout<<"distance and cutoff:"<<distance<<" "<<distcut<<"\n";
            if (AtomCoordinates::dist2(indca[i1],indca[i2]) <distcut2) return distance;
            else return -1.0;
            
        }

    Proximity::Proximity() :natCut(4.5), min_links(2) {}

    Proximity::~Proximity() {}

    void Proximity::set_cutoff(double xx) {distCutoff(xx);}

    double Proximity::get_cutoff() const {return natCut;}

    int Proximity::init(Population *p)
    {
        natCut2=natCut*natCut;
        natoms.resize(p->NumberOfLigands());
        nhatoms.resize(p->NumberOfLigands());
        hatoms.resize(p->NumberOfLigands());

        if (nhatoms.empty()) return 0;

        for (int i=0;i<p->NumberOfLigands();++i) {
            natoms[i]=p->ligand(i)->NumberOfAtoms();
            nhatoms[i]=p->ligand(i)->numHeavyAtoms();
            hatoms[i].resize(nhatoms[i]);
            int k=0;

            for (int ja=0;ja<natoms[i];++ja) {
                if (p->ligand(i)->atom(ja).Species() !=hydrogen) {
                    hatoms[i][k++]=p->ligand(i)->atom(ja).UniqueId();
                }
            }
        }

        return 1;
    }

    bool Proximity::operator()(const Contact &c) const
    {
        return operator()(c.first,c.second);
    }

    bool Proximity::operator()(int i,int j) const
    {
        if (abs(i-j) <minseqsep) return false;

        int cont=0;

        for (int a=0;a<nhatoms[i];a++) {
            for (int b=0;b<nhatoms[j];b++) {
                if (AtomCoordinates::dist2(hatoms[i][a],hatoms[j][b])<natCut2) {
                    cont++;

                    if (cont >= min_links) return true;
                }
            }

        }

        return false;

        //cont is initialized outside the loops. One can only get here if
        //total number of contacts a->b and b->a added is less than 2
    }

    HBContact::HBContact() : ohbmm(NULL), hbcut(1.03), myhb(false)
    {
        issym=false;
    }

    HBContact::HBContact(HBMM *gohbmm, double hbc) : ohbmm(gohbmm),hbcut(hbc),
            myhb(false) {}

    HBContact::~HBContact()
    {
        if (myhb && ohbmm!=NULL) delete(ohbmm);
    }

    void HBContact::set_cutoff(double xx) {setHBCut(xx);}

    double HBContact::get_cutoff() const {return hbcut;}

    int HBContact::init(Population *p)
    {
        if (p==NULL) return 0;

        if (ohbmm==NULL) {
            ohbmm=new HBMM();
            ohbmm->Connect(p);
            ohbmm->init();
            myhb=true;
        }

        return 1;
    }

    bool HBContact::operator()(const Contact &c) const
    {
        return (fabs(ohbmm->InterLg(c.first,c.second)) >hbcut);
    }

    bool HBContact::operator()(int i1,int i2) const
    {
        return (fabs(ohbmm->InterLg(i1,i2)) >hbcut);
    }

    HBContactChains::HBContactChains() : ohbmm(NULL), hbcut(4.6),
            myhb(false) {}

    HBContactChains::HBContactChains(HBMM *gohbmm, double hbc) :
            ohbmm(gohbmm),hbcut(hbc), myhb(false) {}

    HBContactChains::~HBContactChains()
    {
        if (myhb && ohbmm!=NULL) delete(ohbmm);
    }

    void HBContactChains::set_cutoff(double xx) {setHBCut(xx);}

    double HBContactChains::get_cutoff() const {return hbcut;}

    int HBContactChains::init(Population *p)
    {
        if (p==NULL) return 0;

        if (ohbmm==NULL) {
            ohbmm=new HBMM();
            ohbmm->Connect(p);
            ohbmm->init();
            myhb=true;
        }

        return 1;
    }

    bool HBContactChains::operator()(const Contact &c) const
    {
        return (fabs(ohbmm->InterChain(c.first,c.second)) >hbcut);
    }

    bool HBContactChains::operator()(int i1,int i2) const
    {
        return (fabs(ohbmm->InterChain(i1,i2)) >hbcut);
    }

    HPContact::HPContact() : ohpf(NULL), fraccut(0.1), myhp(false) {}

    HPContact::HPContact(Hydrophobicity *gohpf,double ftc) :
            ohpf(gohpf), fraccut(ftc) {}

    HPContact::~HPContact()
    {
        if (myhp && ohpf!=NULL) delete(ohpf);
    }

    void HPContact::set_cutoff(double xx) {setCut(xx);}

    double HPContact::get_cutoff() const {return fraccut;}

    int HPContact::init(Population *p)
    {
        if (p==NULL) return 0;

        if (ohpf==NULL) {
            ohpf=new Hydrophobicity();
            ohpf->Connect(p);
            ohpf->init();
            myhp=true;
        }

        return 1;
    }

    bool HPContact::operator()(const Contact &c)
    {
        return operator()(c.first,c.second);
    }

    bool HPContact::operator()(const Contact &c) const
    {
        return ((ohpf->hp_contact_frac(
                     c.first,
                     c.second,
                     AtomCoordinates::dist2)) >=fraccut);
    }

    bool HPContact::operator()(int i1,int i2) const
    {
        if (abs(i1-i2) <minseqsep) return false;

        return ((ohpf->hp_contact_frac(
                     i1,
                     i2,
                     AtomCoordinates::dist2) >=fraccut));
    }

}

