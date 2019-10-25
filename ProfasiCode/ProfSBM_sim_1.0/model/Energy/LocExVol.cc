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

#include "LocExVol.hh"
using std::valarray;
using std::deque;
using std::pair;
using std::make_pair;
using std::swap;

namespace prf
{

    LocExVol::LocExVol() : ExVolBase(),Energy()
    {
        Name("LocExVol");
        sig2.resize(NA*NA,0);

        for (int i=0;i<NA;++i) for (int j=0;j<NA;++j)
                sig2[i+j*NA]=(sigsa[i]+sigsa[j])*(sigsa[i]+sigsa[j]);

        asa.resize(NA*NA,0);

        bsa.resize(NA*NA,0);

        for (int i=0;i<NA*NA;++i) {
            asa[i]=-7*pow(sig2[i]/cut2,6.0);
            bsa[i]=6*pow(sig2[i]/cut2,6.0)/cut2;
        }
        grdtyp=2;
    }

    LocExVol::~LocExVol() {}

    double LocExVol::Vexv(int ipair)
    {
        int a;
        double r2,r6;

        if ((r2=AtomCoordinates::s2(lci1[ipair],lci2[ipair]))>cut2)
            return 0;

        r6=sig2[(a=lcid[ipair])]/r2;

        r6*=r6*r6;

        return (r6*r6+asa[a]+bsa[a]*r2);
    }

    double LocExVol::dVexv(int ipair, std::valarray<double> &gx)
    {
        int a;
        double r2,r6;

        Vector3 v1=AtomCoordinates::vec(lci1[ipair]);
        Vector3 v2=AtomCoordinates::vec(lci2[ipair]);
        Vector3 v21=v2-v1;
        if ((r2=(v21.mag2()))>cut2) return 0;

        r6=sig2[(a=lcid[ipair])]/r2;

        r6*=r6*r6;

        double ans=(r6*r6+asa[a]+bsa[a]*r2);
        v21*=((-12.0/r2)*(ans-5*bsa[a]/6.0/r2-asa[a]));
        int iloc=3*lci2[ipair];
        gx[iloc++]+=v21.x();
        gx[iloc++]+=v21.y();
        gx[iloc++]+=v21.z();
        iloc=3*lci1[ipair];
        gx[iloc++]-=v21.x();
        gx[iloc++]-=v21.y();
        gx[iloc++]-=v21.z();
        return ans;
    }

    void LocExVol::init()
    {
        if (initialized) return;
        Logger blog;
        int ich, ipair;
        deque<pair<int,int> > retr,allpairs;
        maxdof=-1;
        totdof=0;

        for (ich=0;ich<p->NumberOfChains();++ich) {
            int chdof=p->Chain(ich)->numBBdof()+p->Chain(ich)->numRTdof();

            if (chdof>maxdof) maxdof=chdof;

            totdof+=chdof;
        }

        blog(15)<<"LocExv: Constructing local pair information\n";

        blog<<"found maxdof = "<<maxdof<<" and totdof = "<<totdof<<"\n";
        int nerrors=0,at1,at2,absdofno;
        doftorng.resize(2*NC()*maxdof,-10);
        ipair=0;
        allpairs.clear();

        for (ich=0;ich<p->NumberOfChains();++ich) {
            for (int ibbdof=0;ibbdof<p->Chain(ich)->numBBdof();++ibbdof) {
                p->Chain(ich)->LocPairsBBdof(ibbdof,retr);

                if (retr.size()>9) {
                    prf::cerr<<"bbdof "<<ibbdof<<" for chain "<<ich<<" returned "
                    <<retr.size()<<" local pairs, which is too many\n";
                    ++nerrors;
                }

                absdofno=ich*maxdof+ibbdof;

                doftorng[2*absdofno]=ipair;

                for (size_t p_id=0;p_id<retr.size();++p_id) {
                    at1=retr[p_id].first;
                    at2=retr[p_id].second;

                    if (at2<at1) swap(at1,at2);

                    allpairs.push_back(make_pair(at1,at2));

                    ipair++;
                }

                doftorng[2*absdofno+1]=ipair;
            }

            for (int irtdof=0;irtdof<p->Chain(ich)->numRTdof();++irtdof) {
                p->Chain(ich)->LocPairsRTdof(irtdof,retr);

                if (retr.size()>9) {
                    prf::cerr<<"rtdof "<<irtdof<<" for chain "<<ich<<" returned "
                    <<retr.size()<<" local pairs, which is too many\n";
                    ++nerrors;
                }

                absdofno=ich*maxdof+p->Chain(ich)->numBBdof()+irtdof;

                doftorng[2*absdofno]=ipair;

                for (size_t p_id=0;p_id<retr.size();++p_id) {
                    at1=retr[p_id].first;
                    at2=retr[p_id].second;

                    if (at2<at1) swap(at1,at2);

                    allpairs.push_back(make_pair(at1,at2));

                    ipair++;
                }

                doftorng[2*absdofno+1]=ipair;
            }
        }

        npair=allpairs.size();

        blog(15)<<Name()<<"> Found "<<npair<<" local pairs using DOF information.\n";
        blog(50)<<"They are ...\n";

        for (int i=0;i<npair;++i) {
            blog<<"("<<allpairs[i].first<<", "<<allpairs[i].second<<") ";

            if (((i+1)%10)==0) blog<<"\n";
        }

        blog<<"\n"<<Name()<<"> End of local pair list.\n";

        lci1.resize(npair,-100);
        lci2.resize(npair,-100);
        lcid.resize(npair,-100);
        Melpsa.resize(npair,0.0);
        Velpsa.resize(npair,0.0);

        for (int i=0;i<npair;++i) {
            lci1[i]=allpairs[i].first;
            lci2[i]=allpairs[i].second;
            lcid[i]=p->PairType(lci1[i],lci2[i]);
        }

        if (nerrors!=0) {
            prf::cerr <<"LocExVol::Connect : "<<nerrors<<" errors\nexit\n";
            exit(1);
        }
        initialized=true;
    }

    double LocExVol::evaluate()
    {
        vval=delv=0;

        for (int k=0;k<npair;++k) {
            vval+=(Melpsa[k]=Vexv(k));
        }

        upairbegin=upairend=-1;

        return vval*=ksa;
    }

    double LocExVol::gradientXYZ(std::valarray<double> &ans)
    {
        vval=0;
        ans=0;
        for (int k=0;k<npair;++k) {
            vval+=dVexv(k,ans);
        }
        ans*=ksa;
        return vval*=ksa;
    }


    void LocExVol::rangeEstimate(double &x1, double &x2)
    {
        x1=0.085*NPairs();
        x2=0.12*NPairs();
        //mere guess
    }

    double LocExVol::deltaE(Update *updt)
    {
        delv=0;
        upairbegin=upairend=0;

        if (updt->rigid_chain_update()) return delv;

        for (unsigned i=0;i<updt->num_changes();++i) {
            DOF_Info tmp=updt->change(i).info;
            int ncur=tmp.chain;

            upairbegin=doftorng[2*(ncur*maxdof+tmp.index_in_chain-9)];
            upairend=doftorng[2*(ncur*maxdof+tmp.index_in_chain-9)+1];

            for (int k=upairbegin;k<upairend;++k) {
                delv+=((Velpsa[k]=Vexv(k))-Melpsa[k]);
            }
        }

        return delv*=ksa;
    }

    void LocExVol::Accept(Update *updt)
    {
        if (updt->rigid_chain_update()) return;
        for (unsigned i=0;i<updt->num_changes();++i) {
            DOF_Info tmp=updt->change(i).info;
            int ncur=tmp.chain;

            upairbegin=doftorng[2*(ncur*maxdof+tmp.index_in_chain-9)];
            upairend=doftorng[2*(ncur*maxdof+tmp.index_in_chain-9)+1];

            for (int k=upairbegin;k<upairend;++k) {
                Melpsa[k]=Velpsa[k];
            }
        }

        if (vval>1000 && delv<-0.3*vval) {
            vval=ksa*Melpsa.sum();
        } else vval+=delv;

        upairbegin=upairend=-1;
    }

    void LocExVol::PrintPairs() const
    {
        for (int i=0;i<npair;++i) {
            prf::clog<<i<<" : ("<<lci1[i]<<", "<<lci2[i]<<")\n";
        }
    }

    bool LocExVol::is_loc_pair(int i1, int i2)
    {
        bool ans=false;

        for (int i=0;i<npair;++i) {
            if ((lci1[i]==i1&&lci2[i]==i2) || (lci1[i]==i2&&lci2[i]==i1)) {
                ans=true;
                break;
            }
        }

        return ans;
    }
}
