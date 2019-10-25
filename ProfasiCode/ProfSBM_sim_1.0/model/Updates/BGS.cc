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

#include "BGS.hh"

using namespace UnivConstants;

namespace prf
{
    BGS::BGS() : Update()
    {
        Name("BGS");
        islocal=isbbu=true;
        isscu=isrgu=ismcu=false;
        ranges.resize(5);
        atom_ranges.resize(8);
        nranges=5;
        natomranges=8;
        abgs=300;
        bbgs=10;
        nph=0;
        bgswt=0;
        dirn=0;
        nchanges=8;
        thechange.resize(8);
    }

    BGS::~BGS() {}

    void BGS::build_dof_list()
    {
        std::deque<DOF_Info> mydofs;
        for (std::vector<DOF_Info>::iterator i=popl->dof_map().begin();
        i!=popl->dof_map().end();++i) {
            if (i->dof_kind==backbone_torsion_angle) mydofs.push_back(*i);
        }
        dof_at_site.resize(mydofs.size());
        dof_at_site.assign(mydofs.begin(),mydofs.end());
        site_weight.resize(dof_at_site.size(),1);
        for (size_t i=0;i<dof_at_site.size();++i) {
            DOF_Info ang=dof_at_site[i];
            int chid=ang.chain;
            int room=popl->Chain(chid)->last_AA()->UniqueId()-ang.group;
            if (room<3) site_weight[i]=0;
            else if (ang.specific_index_in_group==1) {
                if (popl->ligand(ang.group)->OLC()==P or
                    popl->ligand(ang.group)->OLC()==DPR) site_weight[i]=2;
                else site_weight[i]=0;
            }
        }
    }

    void BGS::print_setup(std::string &st) {}

    int BGS::perform()
    {
        Update::perform();
        DOF_Info startangle=dof_at_site[site];
        int chid=startangle.chain;

        for (int i=0;i<8;++i) {iph[i]=-50000;dph[i]=0;}

        int i,j,k,ia,upsch;

        Vector3 iv[8],bv[8],rv[3],dv[3][8];
        double ab[8];
        double A[8][8],p[8],psi[8],sum,r1,r2;
        double wfw,wbw;

        int N=popl->Chain(chid)->numAminoAcids();
        ia=startangle.group-popl->Chain(chid)->first_AA()->UniqueId();

        if (ia>(N-4)/2) dirn=0;
        else dirn=1;

        Ligand *curlg=popl->Chain(chid)->AA(ia);

        upsch=curlg->UniqueId();

        if (ia==0) st_fl=popl->Chain(chid)->begin_atom();
        else st_fl=popl->Chain(chid)->AA(ia)->first_atom().UniqueId();

        if (ia==N-4) nd_fl=popl->Chain(chid)->end_atom();
        else nd_fl=popl->Chain(chid)->AA(ia+3)->last_atom().UniqueId()+1;

        if (dirn==0) {
            st_atom=st_fl;
            nd_atom=popl->Chain(chid)->end_atom();

            if (upsch+3>=popl->Chain(chid)->last_ligand()->UniqueId())
                nranges=4;
            else {
                ranges[4].first=upsch+4;
                ranges[4].second=popl->Chain(chid)->last_ligand()->UniqueId();
                nranges=5;
            }

            for (int irng=0;irng<4;++irng)
                ranges[irng].first=ranges[irng].second=upsch+irng;
        } else {
            st_atom=popl->Chain(chid)->begin_atom();
            nd_atom=nd_fl;

            if (upsch==popl->Chain(chid)->first_ligand()->UniqueId())
                nranges=4;
            else {
                ranges[0].first=popl->Chain(chid)->first_ligand()->UniqueId();
                ranges[0].second=upsch-1;
                nranges=5;
            }

            for (int irng=0;irng<4;++irng)
                ranges[nranges-4+irng].first=
                    ranges[nranges-4+irng].second=upsch+irng;
        }

        nph=0;

        for (i=0;i<4;++i) {
            if (popl->Chain(chid)->AA(ia+i)->OLC()!=P &&
                popl->Chain(chid)->AA(ia+i)->OLC()!=DPR) {
                iv[nph]=popl->Chain(chid)->AA(ia+i)->Calpha().Position();
                bv[nph]=iv[nph]-popl->Chain(chid)
                        ->AA(ia+i)->Nitrogen().Position();
                ab[nph]=AminoAcid::b[1];
                iph[nph++]=3*(ia+i);
            }

            iv[nph]=popl->Chain(chid)->AA(ia+i)->Cprime().Position();

            bv[nph]=iv[nph]-popl->Chain(chid)->AA(ia+i)->Calpha().Position();
            ab[nph]=AminoAcid::b[2];
            iph[nph++]=3*(ia+i)+1;
        }

        nchanges=nph;

        rv[0]=popl->Chain(chid)->AA(ia+3)->Calpha().Position();
        rv[1]=popl->Chain(chid)->AA(ia+3)->Cprime().Position();
        rv[2]=popl->Chain(chid)->AA(ia+3)->Oc().Position();

        for (i=0;i<3;++i) {
            for (j=0;j<nph;++j) {
                dv[i][j]=(1.0/ab[j])*(bv[j].cross(rv[i]-iv[j]));
            }
        }

        for (i=0;i<nph;++i) {
            for (j=i;j<nph;++j) {
                A[i][j]=0;

                for (k=0;k<3;++k) {
                    A[i][j]+=dv[k][i].dot(dv[k][j]);
                }

                A[i][j]*=bbgs;

                if (i==j) A[i][j]+=1;

                A[i][j]*=abgs/2;
            }
        }

        for (i=0;i<nph;++i) {
            for (j=i;j<nph;++j) {
                for (sum=A[i][j],k=i-1;k>=0;--k) sum-=A[i][k]*A[j][k];

                if (i==j) p[i]=sqrt(sum);
                else A[j][i]=sum/p[i];
            }
        }

        for (i=0;i<8;i+=2) {
            r1=sqrt(-log(rnd->shoot()));
            r2=rnd->shoot();
            psi[i]=r1*cos(pi2*r2);
            psi[i+1]=r1*sin(pi2*r2);
        }

        for (i=0;i<nph;++i) dph[i]=0;

        for (i=nph-1;i>=0;--i) {
            for (sum=psi[i],k=i+1;k<nph;k++) sum-=A[k][i]*dph[k];

            dph[i]=sum/p[i];
        }

        sum=0;

        for (i=0;i<nph;++i) sum+=psi[i]*psi[i];

        wfw=exp(-sum);

        for (i=0;i<nph;++i) wfw*=p[i];


        for (i=0;i<nph;++i) {
            thechange[i].before=popl->Chain(chid)->backbone()->torsional_angle(iph[i]);
            thechange[i].after=thechange[i].before+dph[i];
            thechange[i].info=dof_at_site[site+i];
            //BGS works on adjacent backbone DOFs
            popl->Chain(chid)->backbone()->incr_tors_angle(iph[i],dph[i]);
        }
        popl->Chain(chid)->reconstruct(dirn,dirn?ia+3:ia);

        /** backward probability **/
        nph=0;

        for (i=0;i<4;++i) {
            if (popl->Chain(chid)->AA(ia+i)->OLC()!=P&&
                popl->Chain(chid)->AA(ia+i)->OLC()!=DPR) {
                iv[nph]=popl->Chain(chid)->AA(ia+i)->Calpha().Position();
                bv[nph]=iv[nph]-popl->Chain(chid)
                        ->AA(ia+i)->Nitrogen().Position();
                ab[nph]=AminoAcid::b[1];
                iph[nph++]=3*(ia+i);
            }

            iv[nph]=popl->Chain(chid)->AA(ia+i)->Cprime().Position();

            bv[nph]=iv[nph]-popl->Chain(chid)->AA(ia+i)->Calpha().Position();
            ab[nph]=AminoAcid::b[2];
            iph[nph++]=3*(ia+i)+1;
        }

        rv[0]=popl->Chain(chid)->AA(ia+3)->Calpha().Position();

        rv[1]=popl->Chain(chid)->AA(ia+3)->Cprime().Position();
        rv[2]=popl->Chain(chid)->AA(ia+3)->Oc().Position();

        for (i=0;i<3;++i) {
            for (j=0;j<nph;++j) {
                dv[i][j]=(1.0/ab[j])*(bv[j].cross(rv[i]-iv[j]));
            }
        }

        for (i=0;i<nph;++i) {
            for (j=i;j<nph;++j) {
                A[i][j]=0;

                for (k=0;k<3;++k) {
                    A[i][j]+=dv[k][i].dot(dv[k][j]);
                }

                A[i][j]*=bbgs;

                if (i==j) A[i][j]+=1;

                A[i][j]*=abgs/2;
            }
        }

        for (i=0;i<nph;++i) {
            for (j=i;j<nph;++j) {
                for (sum=A[i][j],k=i-1;k>=0;--k) sum-=A[i][k]*A[j][k];

                if (i==j) p[i]=sqrt(sum);
                else A[j][i]=sum/p[i];
            }
        }

        for (i=0;i<nph;++i) {
            psi[i]=p[i]*dph[i];

            for (j=i+1;j<nph;++j) psi[i]+=A[j][i]*dph[j];
        }

        sum=0;

        for (i=0;i<nph;++i) sum+=psi[i]*psi[i];

        wbw=exp(-sum);

        for (i=0;i<nph;++i) wbw*=p[i];

        bgswt=wbw/wfw;;

        return 0;
    }

    double BGS::intrinsic_weight() const {return bgswt;}

    int BGS::revert()
    {
        if (Update::revert()) {
            int upch=thechange[0].info.chain;
            int upsch=thechange[0].info.group-popl->Chain(upch)->first_AA()->UniqueId();

            if (dirn == 0) {
                popl->Chain(upch)->backbone()->freconst_bond_vectors(3*upsch);
            } else {
                popl->Chain(upch)->backbone()->breconst_bond_vectors(3*(upsch+3)+2);
            }

            return 1;
        } else
            return 0;
    }
}

