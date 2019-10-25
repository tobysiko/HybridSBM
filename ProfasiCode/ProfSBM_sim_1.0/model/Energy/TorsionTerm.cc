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

#include "TorsionTerm.hh"
#include <cmath>

using namespace std;
using UnivConstants::pi;

namespace prf
{
    TorsionTerm::TorsionTerm() : Energy()
    {
        Name("TorsionTerm");
        kgly=-0.15;
        ntors=0;
        grdtyp=3;
    }

    TorsionTerm::~TorsionTerm() {}

    void TorsionTerm::init()
    {
        if (initialized) return;
        ntors=0;
        for (int i=0;i<p->NumberOfChains();++i)
            ntors+=p->Chain(i)->numBBdof()+p->Chain(i)->numRTdof();

        location.resize(p->n_dof(),-1);
        torstype.resize(ntors,-1);
        dofid.resize(ntors,-1);

        int itors=0;
        for (int i=0;i<p->n_dof();++i) {
            DOF_Info dof=p->get_dof_info(i);
            if (dof.dof_kind==sidechain_torsion_angle) {
                Ligand* res = p->ligand(dof.group);
                OneLetterCode o=res->OLC();
                if ((!res->isAA())||o==A||o==P||o==DPR) continue;
                int locindx=dof.specific_index_in_group;
                switch(locindx) {
                case 0:
                    regtors(i,itors,0); break;
                case 1:
                    if (o==I||o==L||o==M||
                        o==E||o==Q||o==K||o==R) regtors(i,itors,0);
                    if (o==D||o==N) regtors(i,itors,3);
                    if (o==F||o==Y||o==W||o==H) regtors(i,itors,2);
                    break;
                case 2:
                    if (o==M) regtors(i,itors,1);
                    if (o==E || o==Q) regtors(i,itors,3); // CHANGED or TO ||
                    if (o==K || o==R) regtors(i,itors,0); // CHANGED or TO ||
                    break;
                case 3:
                    if (o==K) regtors(i,itors,0);
                    if (o==R) regtors(i,itors,2);
                    break;
                default:
                    break;
                };
            } else if (dof.dof_kind==backbone_torsion_angle ) {
                AminoAcid *res=(AminoAcid *)(p->ligand(dof.group));
                if (res->OLC()==G &&
                    dof.specific_index_in_group==1) {
                    if (!(res->hasNTerminal()||res->hasCTerminal())) {
                        regtors(i,itors,4);
                    }
                }
            }
        }
        std::valarray<int> tmpdofid=dofid[slice(0,itors,1)];
        std::valarray<int> tmptorstype=torstype[slice(0,itors,1)];
        ntors=itors;
        dofid.resize(ntors);
        torstype.resize(ntors);
        dofid=tmpdofid;
        torstype=tmptorstype;
        initialized=true;
    }

    void TorsionTerm::regtors(int i, int &itors,int ityp)
    {
        location[i]=itors;
        dofid[itors]=i;
        torstype[itors]=ityp;
        itors++;
    }

    double TorsionTerm::evaluate()
    {
        delv=vval=0;
        for (int i=0;i<ntors;++i) {
            vval+=torsion(i);
        }
        return vval;
    }

    double TorsionTerm::gradientDOF(std::valarray<double> &ans,
                                 std::vector<int> &indxs)
    {
        double sum=0;
        ans=0;
        for (size_t i=0;i<indxs.size();++i) {
            int sno=location[indxs[i]];
            if (sno>=0 and sno<ntors) {
                ans[i]=dtorsion(sno);
                sum+=torsion(sno);
            }
        }
        return sum;
    }

    double TorsionTerm::deltaE(Update *updt)
    {
        delv=0;
        if (updt->rigid_chain_update()) return 0;
        for (size_t i=0;i<updt->num_changes();++i) {
            dof_change_type c=updt->change(i);
            int sno=location[c.info.global_index];
            if (sno>=0 and sno<ntors) {
                int itype=torstype[sno];
                delv+=(torFcn(c.after,itype)-
                       torFcn(c.before,itype));
            }
        }
        return delv;
    }

    double TorsionTerm::torsion(int sno)
    {
        if (sno<0 or sno>=ntors) return 0;
        int itype=torstype[sno];
        double vl=p->get_dof(dofid[sno]);
        return torFcn(vl,itype);
    }

    double TorsionTerm::dtorsion(int sno)
    {
        if (sno<0 or sno>=ntors) return 0;
        int itype=torstype[sno];
        double vl=p->get_dof(dofid[sno]);
        return dTorFcn(vl,itype);
    }

    void TorsionTerm::rangeEstimate(double &x1, double &x2)
    {
        int npro=0,nala=0,ngly=0,nrtdof=0;

        for (int j=0;j<NC();++j) {
            for (int i=0;i<p->Chain(j)->numAminoAcids();++i) {
                prf::OneLetterCode olc=p->Chain(j)->AA(i)->OLC();
                if (olc==P || olc==DPR) ++npro;
                else if (olc==G) ++ngly;
                else if (olc==A) ++nala;
            }
            nrtdof+=p->Chain(j)->numRTdof();
        }

        x1=-0.4*(nrtdof-nala)+ngly*kgly;
        x2=0.5*x1;
    }

    //there are four different itypes,
    //0: sp3-sp3,
    //1: C-S in Met,
    //2: sp3-sp2 and sp3-aromatic,
    //3: sp3-amide and sp3-CO2 (D,E,N,Q)
    double TorsionTerm::torFcn(double phi, int itype)
    {
        const double k[4]={0.6,0.3,0.4,-0.4};
        const double cst[4] = {3, 3, 2, 2};
        if (itype==4) return kgly*(cos(phi)+2*cos(2*phi));
        return k[itype]*cos(cst[itype]*phi);
    }

    double TorsionTerm::dTorFcn(double phi, int itype)
    {
        const double cstk[4]={0.6*3,0.3*3,0.4*2,-0.4*2};
        const double cst[4] = {3, 3, 2, 2};
        if (itype==4) return -kgly*(sin(phi)+4*sin(2*phi));
        return -cstk[itype]*sin(cst[itype]*phi);
    }
}

