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

#include "ExVol.hh"

namespace prf
{
    ExVol::ExVol() : ExVolBase(), Energy()
    {
        Name("ExVol");
        LAMBDA=0.75;
        cutg=LAMBDA*cut;
        cutg2=LAMBDA*LAMBDA*cut*cut;
        sig2.resize(NA*NA,0);

        for (int i=0;i<NA;++i) for (int j=0;j<NA;++j)
                sig2[i+j*NA]=sqr(LAMBDA*(sigsa[i]+sigsa[j]));

        asa.resize(NA*NA,0);

        bsa.resize(NA*NA,0);

        for (int i=0;i<NA*NA;++i) {
            asa[i]=-7*pow(sig2[i]/cutg2,6.0);
            bsa[i]=6*pow(sig2[i]/cutg2,6.0)/cutg2;
        }

        va_ngbdisp.resize(27,0);
        enew=0;
        grdtyp=2;
    }

    ExVol::~ExVol() {}

    void ExVol::init()
    {
        if (initialized) return;
        Logger blog;

        xmax=0.5*AtomCoordinates::boxL();
        xmin=-xmax;
        neginf=-5.5*AtomCoordinates::boxL();
        ymin=zmin=xmin;
        ymax=zmax=xmax;
        nx=(int)((xmax-xmin)/cutg);
        ny=(int)((ymax-ymin)/cutg);
        nz=(int)((zmax-zmin)/cutg);
        bxl=AtomCoordinates::boxL();
        cellx=bxl/((double)nx);
        celly=bxl/((double)ny);
        cellz=bxl/((double)nz);
        hbxl=0.5*bxl;
        cellxc=bxl-cellx;
        cellyc=bxl-celly;
        cellzc=bxl-cellz;

        xpitch=1;
        ypitch=nx;
        zpitch=nx*ny;
        /* how many cell indices one moves by to reach the next cell along x,y,z */
        nc=nx*ny*nz;
        cell.resize(nc,0);
        cellstatus=-1;
        /***********************************/
        /* The part inside the above two C style comment lines used to be
          in the constructor. But when used with an interface class, the
          constructor is called along with the Interface class constructor,
          where as the box length is only set when the interface class is
          constructed and then the parameters reassigned or read from a file.
          Using the above lines in the constructor therefore, will only use
          the default values of box length etc which is not such a good idea. */

        NTO=0;

        for (int i=0;i<p->NumberOfChains();++i)
            NTO+=p->Chain(i)->numberOfAtoms();

        MAXNGB=std::min(500,NTO);

        listt.resize(NTO,0);

        pnt.resize(NTO,0);;

        lxyz.resize(3*NTO,0.0);

        blxyz.resize(3*NTO,0.0);

        atomloc.resize(NTO,0);

        atlcbk.resize(NTO,0);

        atom_stat_bkp.resize(NTO,0);

        atom_stat.resize(NTO,0);

        atom_used.resize(NTO,0);

        unmoved_ngb.resize(MAXNGB,0);

        rigid_ngb.resize(MAXNGB,0);

        flxbl_ngb.resize(MAXNGB,0);

        connected.NAtoms(NTO);

        for (int i=0;i<NC();++i) {
            blog(15)<<"ExVol: Importing connections information from peptide "
            <<i<<" : "<<(p->PepName(i)).c_str()<<"\n";
            p->Chain(i)->ExportConnections(connected);
        }

        blog<<"Number of cells = "<<nc<<"\n";

        blog<<"cell linear dimension = "<<cellx<<"\n";
        blog<<"cell pitch = "<<xpitch<<", "<<ypitch<<", "<<zpitch<<"\n";
        blog<<"maximum number of neighbours allowed = "<<MAXNGB<<"\n";
        initialized=true;
    }

    void ExVol::rangeEstimate(double &x1, double &x2)
    {
        x1=0;
        x2=0.03*AtomCoordinates::numberOfAtoms();
    }

    double ExVol::evaluate()
    {
        initcells();
        enew=vval=ksa*calc_esa(1e50);
        imprequest=0;
        //cout <<"excluded volume: calc_esa returns "<<vval<<"\n";
        //prf::cout <<"sac full returns " <<sacfull()<<"\n";
        //getchar();
        delv=0;
        return vval;
    }

    double ExVol::gradientXYZ(std::valarray<double> &ans)
    {
        initcells();
        ans=0;
        vval=ksa*calc_esa_with_grd(ans);
        ans*=ksa;
        return vval;
    }

    double ExVol::deltaE(Update *updt)
    {
        return deltaEwithlimit(updt,1e50);
    }

    double ExVol::deltaEwithlimit(Update *updt, double emax)
    {
        double etest,etest1;
        double epart1,epart2;
//     epart1=vval;
//     epart2=evaluate();
//     vval=epart1;delv=epart2-epart1;
//     return delv;

        if (emax<-vval) {imprequest=1;return 0;}
        else imprequest=0;

        etest=emax/ksa;

        //cout<<"update involves "<<(updt->end_atom()-updt->begin_atom())
        //<<" atoms out of "<<NTO<<"\n";
        delv=0;

        if (updt->end_atom()-updt->begin_atom()<(NTO/2)) {
            //cout<<"partial esa method\n";
            //getchar();
            mark_moved_atoms(updt->begin_atom(),updt->end_atom(),
                             updt->begin_flexible_part(),
                             updt->end_flexible_part());
            AtomCoordinates::switch_to_backups();
            epart1=ksa*partial_esa(updt->begin_atom(),updt->end_atom(),1e50);
            etest1=emax+epart1;
            AtomCoordinates::switch_to_regular();
            updatecells(updt->begin_atom(),updt->end_atom());
            epart2=ksa*partial_esa(updt->begin_atom(),updt->end_atom(),
                                   etest1/ksa);
            //epart2=ksa*partial_esa(updt->begin_atom(),updt->end_atom(),1e50);
            delv=epart2-epart1;
            enew=vval+delv;

            if (vval>100 && (delv<(-0.3*vval))) {
//                prf::cout<<"vval = "<<vval<<" delv = "
//                <<delv<<" falling back to calc_esa\n";
                enew=ksa*calc_esa(vval/ksa+etest);
                delv=enew-vval;
//                prf::cout<<"got enew = "<<enew<<", delv "<<delv<<", eold = "<<vval<<"\n";
            }

            //checksubsets(updt);
            //getchar();
        } else {
            //cout<<"calc_esa \n";
            //getchar();
            updatecells(updt->begin_atom(),updt->end_atom());
            enew=ksa*calc_esa(vval/ksa+etest);
//            prf::cout<<"Direct calc esa method : enew = "<<enew<<", eold = "<<vval<<", delv = "<<delv<<"\n";
            delv=enew-vval;
        }
        return delv;
    }

    void ExVol::Accept(Update *updt)    //dont accept impossible requests
    {
        if (imprequest!=0) {
            prf::cerr<<"error: accept called when imprequest is "
            <<imprequest<<"\n";
            exit(1);
        }

        vval=enew;

        cellstatus=0;

        if (updt->end_atom()-updt->begin_atom()<(NTO/2))
            reset_markings(updt->begin_atom(), updt->end_atom());
    }

    void ExVol::Revert(Update * updt)
    {
        if (imprequest==1) return;

        if (updt->end_atom()-updt->begin_atom()<(NTO/2))
            reset_markings(updt->begin_atom(), updt->end_atom());
        enew=vval;
        restorecells(updt->begin_atom(), updt->end_atom());
    }

    void ExVol::mark_moved_atoms(int i, int j, int k, int l)
    {
        int ii;

        for (ii=i;ii<j;++ii) atom_stat_bkp[ii]=atom_stat[ii]=1;  //rigid subset

        /*   printf("marked atoms %d to %d as rigid movers \n",i,j); */
        for (ii=k;ii<l;++ii) atom_stat_bkp[ii]=atom_stat[ii]=2;  //flexible subset

        /*   printf("marked atoms %d to %d as nonrigid movers \n",k,l); */
        //getchar();
    }

    void ExVol::restore_markings(int i, int j)
    {
        int ii;

        for (ii=i;ii<j;++ii) atom_stat[ii]=atom_stat_bkp[ii];
    }

    void ExVol::reset_markings(int i, int j)
    {
        int ii;

        for (ii=i;ii<j;++ii) atom_stat[ii]=atom_stat_bkp[ii]=0;
    }

    int ExVol::initcells()
    {
        int i,ic;
        //cout <<"found cell of size "<<(cell.size())<<" in initcells\n";

        if (cellstatus==-1) {
            for (i=0;i<nc;++i) cell[i]=-1;   //mark all cells empty
        } else {
            for (i=0;i<NTO;++i) {
                cell[atomloc[i]]=-1;
            }
        }

        for (i=0;i<NTO;++i) {
            pnt[i]=-2000;
            atomloc[i]=-2000;
            atlcbk[i]=-2000;
            atom_stat[i]=atom_stat_bkp[i]=0;
        }

        /*pnt=nothing,  atom loc=nowhere atom status = neutral */
        for (i=0;i<NTO;++i) {   //initial cell assignment to all atoms
            int ii=3*i;
            lxyz[ii]=boxcrdx(i);
            lxyz[ii+1]=boxcrdy(i);
            lxyz[ii+2]=boxcrdz(i);
            ic=((int)(lxyz[ii]/cellx))+nx*(((int)(lxyz[ii+1]/celly))+
                                           ny*((int)(lxyz[ii+2]/cellz)));

            if (ic<0) {
                prf::cerr<<"init cells: calculated negative cell index "
                <<ic<<" for atom "<<i<<"\n";
                //getchar();
            }

            pnt[i]=cell[ic];

            cell[ic]=i;
            atomloc[i]=ic;
        }

        cellstatus=0;//initialized with atoms

        return 1;
    }

    int ExVol::updatecells(int istr, int iend)
    {
        int i,j,ic,icprev,cell_occ;
        //printf("update cells for atoms between %d and %d \n",istr,iend);

        for (i=istr;i<iend;++i) {   //Move stuff around in cells
            int ii=3*i;
            blxyz[ii]=lxyz[ii];
            blxyz[ii+1]=lxyz[ii+1];
            blxyz[ii+2]=lxyz[ii+2];
            lxyz[ii]=boxcrdx(i);
            lxyz[ii+1]=boxcrdy(i);
            lxyz[ii+2]=boxcrdz(i);
            ic=((int)(lxyz[ii]/cellx))+nx*(((int)(lxyz[ii+1]/celly))+
                                           ny*((int)(lxyz[ii+2]/cellz)));
            atlcbk[i]=icprev=atomloc[i];

            if (ic<0) {
                prf::cerr<<"update cells: calculated negative cell index "
                <<ic<<" for atom "<<i<<"\n"
                <<"coordinates :"<<AtomCoordinates::vec(i)<<"\n"
                <<"(ix,iy,iz)= ("<<lxyz[ii]<<", "<<lxyz[ii+1]<<", "<<lxyz[ii+2]<<"\n";
                //getchar();
            }

            cell_occ=cell[icprev];

            if (cell_occ==i) {
                cell[icprev]=pnt[i];
            } else {
                while ((j=pnt[cell_occ])!=i && j>=0) cell_occ=j;

                if (j>=0) pnt[cell_occ]=pnt[i];
                else {
                    prf::cerr<<"strange situation. While updating cells for atom "
                    <<i<<" found new cell index "<<ic
                    <<"but even though atomloc was "<<atomloc[i]
                    <<" the atom "<<i<<"was not found in cell "<<atomloc[i]
                    <<"contents of cell "<<atomloc[i]<<" ...\n";
                    j=cell[atomloc[i]];

                    do {
                        prf::cerr <<j<<"\n";
                        j=pnt[j];
                    } while (j>=0);

                    exit(1);
                }
            }

            pnt[i]=cell[ic];

            cell[ic]=i;
            atomloc[i]=ic;
        }

        cellstatus=1;

        return 1;
    }

    int ExVol::restorecells(int istr, int iend)
    {
        int i,ic,icprev,cell_occ;

        if (cellstatus==1) {
            for (i=istr;i<iend;++i) {   //Move stuff around in cells
                int ii=3*i;
                ic=atlcbk[i];
                icprev=atomloc[i];

                if (ic<0) {
                    prf::cerr<<"restore cells: calculated negative cell index "
                    <<ic<<" for atom "<<i<<"\n";
                    //getchar();
                }

                lxyz[ii]=blxyz[ii];

                lxyz[ii+1]=blxyz[ii+1];
                lxyz[ii+2]=blxyz[ii+2];
                ic=atlcbk[i];
                cell_occ=cell[icprev];

                if (cell_occ==i) {
                    cell[icprev]=pnt[i];
                } else {
                    while (pnt[cell_occ]!=i) cell_occ=pnt[cell_occ];

                    pnt[cell_occ]=pnt[i];
                }

                pnt[i]=cell[ic];

                cell[ic]=i;
                atomloc[i]=ic;
            }

            cellstatus=0;

            return 1;
        }

        return 0;
    }

    /**********************************************************************/
    double ExVol::cellcalc(int inc,int nl,int &incnr)
    {
        int i,j,n,m,a,ncount=0;
        double e=0;
        double r2,r6;
        //double r2a,r2b;

        if (useperiod==0) {
            for (n=0;n<inc;++n) {
                i=listt[n];

                for (m=n+1;m<nl;++m) {
                    j=listt[m];

                    if ((r2=dist2a(i,j))>cutg2) continue;

                    if ((a=PairType(i,j))<0) continue;

                    ++ncount;

                    r6=sig2[a]/r2;

                    r6=r6*r6*r6;

                    e+=r6*r6+asa[a]+bsa[a]*r2;

//    prf::clog<<"Contributing pair "<<i<<", "<<j<<": "
//       <<r6*r6+asa[a]+bsa[a]*r2<<"\n";
                }
            }
        } else {
            for (n=0;n<inc;++n) {
                i=listt[n];

                for (m=n+1;m<nl;++m) {
                    j=listt[m];

                    if ((r2=dist2b(i,j))>cutg2) continue;

                    if ((a=PairType(i,j))<0) continue;

                    ++ncount;

                    r6=sig2[a]/r2;

                    r6=r6*r6*r6;

                    e+=r6*r6+asa[a]+bsa[a]*r2;

//    prf::clog<<"Contributing pair "<<i<<", "<<j<<": "
//       <<r6*r6+asa[a]+bsa[a]*r2<<"\n";
                }
            }
        }

        incnr=ncount;

        return e;
    }

    double ExVol::cellcalc_with_grd(int inc,int nl,int &incnr,
                                    std::valarray<double> &gx)
    {
        int i,j,n,m,a,ncount=0;
        double e=0,term=0;
        double r2,r6;
        Vector3 v21;

        if (useperiod==0) {
            for (n=0;n<inc;++n) {
                i=listt[n];
                for (m=n+1;m<nl;++m) {
                    j=listt[m];
                    v21=vecitoj_a(i,j);
                    if ((r2=(v21.mag2()))>cutg2) continue;
                    if ((a=PairType(i,j))<0) continue;

                    r6=sig2[a]/r2;
                    r6*=r6*r6;

                    term=(r6*r6+asa[a]+bsa[a]*r2);
                    v21*=((-12.0/r2)*(term-5*bsa[a]/6.0/r2-asa[a]));
                    int iloc=3*j;
                    gx[iloc++]+=v21.x();
                    gx[iloc++]+=v21.y();
                    gx[iloc++]+=v21.z();
                    iloc=3*i;
                    gx[iloc++]-=v21.x();
                    gx[iloc++]-=v21.y();
                    gx[iloc++]-=v21.z();

                    ++ncount;
                    e+=term;
                }
            }
        } else {
            for (n=0;n<inc;++n) {
                i=listt[n];
                for (m=n+1;m<nl;++m) {
                    j=listt[m];
                    v21=vecitoj_b(i,j);
                    if ((r2=(v21.mag2()))>cutg2) continue;
                    if ((a=PairType(i,j))<0) continue;

                    r6=sig2[a]/r2;
                    r6*=r6*r6;

                    term=(r6*r6+asa[a]+bsa[a]*r2);
                    v21*=((-12.0/r2)*(term-5*bsa[a]/6.0/r2-asa[a]));
                    int iloc=3*j;
                    gx[iloc++]+=v21.x();
                    gx[iloc++]+=v21.y();
                    gx[iloc++]+=v21.z();
                    iloc=3*i;
                    gx[iloc++]-=v21.x();
                    gx[iloc++]-=v21.y();
                    gx[iloc++]-=v21.z();

                    ++ncount;
                    e+=term;
                }
            }
        }

        incnr=ncount;

        return e;
    }
    double ExVol::sacfull()
    {
        int i,j,a,ncontr=0,nrr=0,nru=0,nuu=0,nrf=0,nuf=0,nff=0;
        double e=0,er=0,eu=0,eru=0,de=0,erf=0,euf=0,eff=0;
        double r2,r6;

        for (i=0;i<NTO;++i) {
            for (j=i+1;j<NTO;++j) {
                if ((r2=AtomCoordinates::dist2(i,j))>cutg2) continue;

                if ((a=PairType(i,j))<0) continue;

                ++ncontr;

                //printf("CONTRIBUTING PAIR %d (%f,%f,%f) and %d (%f,%f,%f): %f\n",i,
                //lxyz[3*i],lxyz[3*i+1],lxyz[3*i+2],j,lxyz[3*j],lxyz[3*j+1],lxyz[3*j+2], r2);
                if (r2<1e-6)
                    prf::cerr<<"ExVol: distance = "<<r2<<" between pair "<<i
                    <<" and "<<j<<"\n";

                r6=sig2[a]/r2;

                r6=r6*r6*r6;

                de=r6*r6+asa[a]+bsa[a]*r2;

                e+=de;

                if (atom_stat[i]==0 && atom_stat[j]==0) {eu+=de;++nuu;}

                if (atom_stat[i]==0 && atom_stat[j]==1) {eru+=de;++nru;}

                if (atom_stat[i]==1 && atom_stat[j]==0) {eru+=de;++nru;}

                if (atom_stat[i]==1 && atom_stat[j]==1) {er+=de;++nrr;}

                if (atom_stat[i]==2 && atom_stat[j]==2) {eff+=de;++nff;}

                if (atom_stat[i]==1 && atom_stat[j]==2) {erf+=de;++nrf;}

                if (atom_stat[i]==2 && atom_stat[j]==1) {erf+=de;++nrf;}

                if (atom_stat[i]==0 && atom_stat[j]==2) {euf+=de;++nuf;}

                if (atom_stat[i]==2 && atom_stat[j]==0) {euf+=de;++nuf;}

            }
        }

        eu*=ksa;

        eru*=ksa;
        er*=ksa;
        e*=ksa;
        erf*=ksa;
        euf*=ksa;
        eff*=ksa;
        printf("sac full diagnosis...(excluding local sa) \n");
        printf("rigid - rigid = %f from %d contributions\n",er,nrr);
        printf("rigid - unmoved = %f from %d contributions\n",eru,nru);
        printf("unmoved - unmoved = %f from %d contributions\n",eu,nuu);
        printf("flex - rigid = %f from %d contributions\n",erf,nrf);
        printf("flex - unmoved = %f from %d contributions\n",euf,nuf);
        printf("flex - flex = %f from %d contributions\n",eff,nff);
        printf("total = %f from %d contributions\n",e,ncontr);
        return e;
    }



    /**********************************************************************/

    double ExVol::partial_esa(int istrt, int iend,double etest)
    {
        int i,j,ii,ai,aj;
        double e=0,eru=0,efu=0,efr=0,eff=0,erf=0;

        for (i=istrt;i<iend;++i) {
            if (atom_stat[i]>2) continue;

            select_neighbours(i);

            /* Evaluate self avoidance energy term for all unmoved neighbours */
            /* with the mobile part in the same cell as atom i */
            if (useperiod==0) {
                for (j=0;j<n_unmoved;++j) {
                    aj=unmoved_ngb[j];

                    for (ii=0;ii<inc_rigid;++ii) {
                        eru+=Vexva(rigid_ngb[ii],aj);

                        //printf("rigid atom %d and unmoved %d : de= %f\n",
                        //ai,aj,de);
                    }

                    for (ii=0;ii<inc_flxbl;++ii) {
                        efu+=Vexva(flxbl_ngb[ii],aj);

                        //printf("flexible atom %d and unmoved %d : de= %f\n",
                        //ai,aj,de);
                    }
                }

                for (j=0;j<n_flxbl;++j) {
                    aj=flxbl_ngb[j];

                    for (ii=0;ii<inc_rigid;++ii) {
                        efr+=Vexva(rigid_ngb[ii], aj);

                        //printf("rigid atom %d and flexible %d : de= %f\n",
                        //ai,aj,de);
                    }
                }

                for (j=inc_rigid;j<n_rigid;++j) {
                    aj=rigid_ngb[j];

                    for (ii=0;ii<inc_flxbl;++ii) {
                        erf+=Vexva(flxbl_ngb[ii],aj);

                        //printf("flex atom %d and rigid %d : de= %f\n",ai,aj,de);
                    }
                }

                for (ii=0;ii<inc_flxbl;++ii) {
                    ai=flxbl_ngb[ii];

                    for (j=ii+1;j<n_flxbl;++j) {
                        eff+=Vexva(ai,flxbl_ngb[j]);

                        //printf("flexible atom %d and flexible %d : de= %f\n",
                        //ai,aj,de);
                    }
                }

                if (e>etest) break;

                for (j=0;j<inc_rigid;++j) atom_stat[rigid_ngb[j]]=3;

                for (j=0;j<inc_flxbl;++j) atom_stat[flxbl_ngb[j]]=4;
            } else {

                for (j=0;j<n_unmoved;++j) {
                    aj=unmoved_ngb[j];

                    for (ii=0;ii<inc_rigid;++ii) {
                        eru+=Vexvb(aj,rigid_ngb[ii]);

                        //printf("rigid atom %d and unmoved %d : de= %f\n",
                        //ai,aj,de);
                    }

                    for (ii=0;ii<inc_flxbl;++ii) {
                        efu+=Vexvb(aj,flxbl_ngb[ii]);

                        //printf("flexible atom %d and unmoved %d : de= %f\n",
                        //ai,aj,de);
                    }
                }

                for (j=0;j<n_flxbl;++j) {
                    aj=flxbl_ngb[j];

                    for (ii=0;ii<inc_rigid;++ii) {
                        efr+=Vexvb(aj,rigid_ngb[ii]);

                        //printf("rigid atom %d and flexible %d : de= %f\n",
                        //ai,aj,de);
                    }
                }

                for (j=inc_rigid;j<n_rigid;++j) {
                    aj=rigid_ngb[j];

                    for (ii=0;ii<inc_flxbl;++ii) {
                        erf+=Vexvb(aj,flxbl_ngb[ii]);

                        //printf("flex atom %d and rigid %d : de= %f\n",ai,aj,de);
                    }
                }

                for (ii=0;ii<inc_flxbl;++ii) {
                    ai=flxbl_ngb[ii];

                    for (j=ii+1;j<n_flxbl;++j) {
                        eff+=Vexvb(ai,flxbl_ngb[j]);

                        //printf("flexible atom %d and flexible %d : de= %f\n",
                        //ai,aj,de);
                    }
                }
            }

            e=eru+efr+erf+eff+efu;

            if (e>etest) break;

            for (j=0;j<inc_rigid;++j) atom_stat[rigid_ngb[j]]=3;

            for (j=0;j<inc_flxbl;++j) atom_stat[flxbl_ngb[j]]=4;
        }

        restore_markings(istrt,iend);

//        e=eru+efr+erf+eff+efu;
        return e;
    }

    double ExVol::calc_esa(double etest)
    {
        int i,j,a;
        int h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13;
        int dxplus,dyplus,dyminus,dzplus,dzminus;
        int in_cell,nl;
        double e=0;
        int ncontr=0,incontr=0;

        for (i=0;i<NTO;++i) atom_used[i]=0;

        i=0;

        while (i<NTO && e<etest) {
            useperiod=0;

            if (atom_used[i]==1) {++i;continue;}

            nl=0;

            listt[nl++]=a=cell[(j=atomloc[i])];
            atom_used[a]=1;

            while ((a=pnt[a])>=0) {listt[nl++]=a;atom_used[a]=1;}

            in_cell=nl;

            /* In the above, an atom was declared used if it was in the same
             cell as the atom i. This will not be done for atoms in neighbouring
            cells */

            /* Find Neighbours :: */
            int ii=3*i;

            if (lxyz[ii+1]<celly) {dyminus=ny-1;useperiod=1;}
            else dyminus=-1;

            if (lxyz[ii+1]>cellyc) {dyplus=1-ny;useperiod=1;}
            else dyplus=1;

            if (lxyz[ii+2]<cellz) {dzminus=nz-1;useperiod=1;}
            else dzminus=-1;

            if (lxyz[ii+2]>cellzc) {dzplus=1-nz;useperiod=1;}
            else dzplus=1;

            if (lxyz[ii]>cellxc) {dxplus=1-nx;useperiod=1;}
            else dxplus=1;

            /* Find Neighbours :: h1: increments 1,0,0 */
            h1=xpitch*dxplus;

            if ((a=cell[j+h1])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h2: increments 0,1,0 */
            h2=ypitch*dyplus;

            if ((a=cell[j+h2])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h3: increments 0,0,1 */
            h3=zpitch*dzplus;

            if ((a=cell[j+h3])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h4: increments 1,1,0 */
            h4=xpitch*dxplus+ypitch*dyplus;

            if ((a=cell[j+h4])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h5: increments 1,-1,0 */
            h5=xpitch*dxplus+ypitch*dyminus;

            if ((a=cell[j+h5])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h6: increments 1,0,1 */
            h6=xpitch*dxplus+zpitch*dzplus;

            if ((a=cell[j+h6])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h7: increments 1,0,-1 */
            h7=xpitch*dxplus+zpitch*dzminus;

            if ((a=cell[j+h7])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h8: increments 0,1,1 */
            h8=ypitch*dyplus+zpitch*dzplus;

            if ((a=cell[j+h8])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h9: increments 0,1,-1 */
            h9=ypitch*dyplus+zpitch*dzminus;

            if ((a=cell[j+h9])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h10: increments 1,1,1 */
            h10=xpitch*dxplus+ypitch*dyplus+zpitch*dzplus;

            if ((a=cell[j+h10])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h11: increments 1,1,-1 */
            h11=xpitch*dxplus+ypitch*dyplus+zpitch*dzminus;

            if ((a=cell[j+h11])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h12: increments 1,-1,1 */
            h12=xpitch*dxplus+ypitch*dyminus+zpitch*dzplus;

            if ((a=cell[j+h12])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h13: increments 1,-1,-1 */
            h13=xpitch*dxplus+ypitch*dyminus+zpitch*dzminus;

            if ((a=cell[j+h13])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Evaluate self avoidance energy term for atoms in the same
            cell as i, with all atoms in the same cell and the neighbouring
            cells. */
            e+=cellcalc(in_cell,nl,incontr);

            ++i;

            ncontr+=incontr;
        }

        //printf("calc_esa e = %f number of contributions = %d\n",e,ncontr);
        return e;
    }

    double ExVol::calc_esa_with_grd(std::valarray<double> &gx)
    {
        int i,j,a;
        int h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13;
        int dxplus,dyplus,dyminus,dzplus,dzminus;
        int in_cell,nl;
        double e=0;
        int ncontr=0,incontr=0;

        atom_used=0;

        i=0;

        while (i<NTO) {
            useperiod=0;

            if (atom_used[i]==1) {++i;continue;}

            nl=0;

            listt[nl++]=a=cell[(j=atomloc[i])];
            atom_used[a]=1;

            while ((a=pnt[a])>=0) {listt[nl++]=a;atom_used[a]=1;}

            in_cell=nl;

            /* In the above, an atom was declared used if it was in the same
             cell as the atom i. This will not be done for atoms in neighbouring
            cells */

            /* Find Neighbours :: */
            int ii=3*i;

            if (lxyz[ii+1]<celly) {dyminus=ny-1;useperiod=1;}
            else dyminus=-1;

            if (lxyz[ii+1]>cellyc) {dyplus=1-ny;useperiod=1;}
            else dyplus=1;

            if (lxyz[ii+2]<cellz) {dzminus=nz-1;useperiod=1;}
            else dzminus=-1;

            if (lxyz[ii+2]>cellzc) {dzplus=1-nz;useperiod=1;}
            else dzplus=1;

            if (lxyz[ii]>cellxc) {dxplus=1-nx;useperiod=1;}
            else dxplus=1;

            /* Find Neighbours :: h1: increments 1,0,0 */
            h1=xpitch*dxplus;

            if ((a=cell[j+h1])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h2: increments 0,1,0 */
            h2=ypitch*dyplus;

            if ((a=cell[j+h2])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h3: increments 0,0,1 */
            h3=zpitch*dzplus;

            if ((a=cell[j+h3])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h4: increments 1,1,0 */
            h4=xpitch*dxplus+ypitch*dyplus;

            if ((a=cell[j+h4])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h5: increments 1,-1,0 */
            h5=xpitch*dxplus+ypitch*dyminus;

            if ((a=cell[j+h5])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h6: increments 1,0,1 */
            h6=xpitch*dxplus+zpitch*dzplus;

            if ((a=cell[j+h6])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h7: increments 1,0,-1 */
            h7=xpitch*dxplus+zpitch*dzminus;

            if ((a=cell[j+h7])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h8: increments 0,1,1 */
            h8=ypitch*dyplus+zpitch*dzplus;

            if ((a=cell[j+h8])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h9: increments 0,1,-1 */
            h9=ypitch*dyplus+zpitch*dzminus;

            if ((a=cell[j+h9])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h10: increments 1,1,1 */
            h10=xpitch*dxplus+ypitch*dyplus+zpitch*dzplus;

            if ((a=cell[j+h10])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h11: increments 1,1,-1 */
            h11=xpitch*dxplus+ypitch*dyplus+zpitch*dzminus;

            if ((a=cell[j+h11])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h12: increments 1,-1,1 */
            h12=xpitch*dxplus+ypitch*dyminus+zpitch*dzplus;

            if ((a=cell[j+h12])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Find Neighbours :: h13: increments 1,-1,-1 */
            h13=xpitch*dxplus+ypitch*dyminus+zpitch*dzminus;

            if ((a=cell[j+h13])>=0) {
                listt[nl++]=a;

                while ((a=pnt[a])>=0) listt[nl++]=a;
            }

            /* Evaluate self avoidance energy term for atoms in the same
            cell as i, with all atoms in the same cell and the neighbouring
            cells. */
            e+=cellcalc_with_grd(in_cell,nl,incontr,gx);

            ++i;

            ncontr+=incontr;
        }

        //printf("calc_esa e = %f number of contributions = %d\n",e,ncontr);
        return e;
    }
    /********************************************************************/
    /********************************************************************/

    void ExVol::select_neighbours(int i)
    {
        int j,a,ngbno;
        int dxplus,dxminus,dyplus,dyminus,dzplus,dzminus;
        inc_unmoved=inc_rigid=inc_flxbl=n_unmoved=n_flxbl=n_rigid=0;
        useperiod=0;
        a=cell[(j=atomloc[i])];
        //cout<<"found atom "<<i<<" in cell "<<j<<"\n";
        //cout <<" which was occupied by atom "<<a<<"\n";

        if (atom_stat[a]==0) unmoved_ngb[n_unmoved++]=a;
        else if (atom_stat[a]==1) rigid_ngb[n_rigid++]=a;
        else if (atom_stat[a]==2) flxbl_ngb[n_flxbl++]=a;

        //cout <<"atom_stat of "<<a<<" was found to be "<<atom_stat[a]<<"\n";
        while ((a=pnt[a])>=0) {
            //cout <<"current atom is "<<a<<"\n";
            if (atom_stat[a]==0) unmoved_ngb[n_unmoved++]=a;
            else if (atom_stat[a]==1) rigid_ngb[n_rigid++]=a;
            else if (atom_stat[a]==2) flxbl_ngb[n_flxbl++]=a;

            //cout <<"atom_stat of "<<a<<" was found to be "<<atom_stat[a]<<"\n";
        }

        inc_unmoved=n_unmoved;

        inc_rigid=n_rigid;
        inc_flxbl=n_flxbl;

        /* Find more neighbours :: */
        int ii=3*i;

        if (lxyz[ii]<cellx) {dxminus=(nx-1)*xpitch;useperiod=1;}
        else dxminus=-xpitch;

        if (lxyz[ii]>cellxc) {dxplus=(1-nx)*xpitch;useperiod=1;}
        else dxplus=xpitch;

        if (lxyz[++ii]<celly) {dyminus=(ny-1)*ypitch;useperiod=1;}
        else dyminus=-ypitch;

        if (lxyz[ii]>cellyc) {dyplus=(1-ny)*ypitch;useperiod=1;}
        else dyplus=ypitch;

        if (lxyz[++ii]<cellz) {dzminus=(nz-1)*zpitch;useperiod=1;}
        else dzminus=-zpitch;

        if (lxyz[ii]>cellzc) {dzplus=(1-nz)*zpitch;useperiod=1;}
        else dzplus=zpitch;

        va_ngbdisp[1]=dxplus ;

        va_ngbdisp[2]=dyplus ;

        va_ngbdisp[3]=dzplus ;

        va_ngbdisp[4]=dxminus ;

        va_ngbdisp[5]=dyminus ;

        va_ngbdisp[6]=dzminus ;

        va_ngbdisp[7]=dxplus+dyplus ;

        va_ngbdisp[8]=dxplus+dyminus ;

        va_ngbdisp[9]=dxminus+dyminus ;

        va_ngbdisp[10]=dxminus+dyplus ;

        va_ngbdisp[11]=dxplus+dzplus ;

        va_ngbdisp[12]=dxplus+dzminus ;

        va_ngbdisp[13]=dxminus+dzminus ;

        va_ngbdisp[14]=dxminus+dzplus ;

        va_ngbdisp[15]=dyplus+dzplus ;

        va_ngbdisp[16]=dyplus+dzminus ;

        va_ngbdisp[17]=dyminus+dzminus ;

        va_ngbdisp[18]=dyminus+dzplus ;

        va_ngbdisp[19]=dxplus+dyplus+dzplus ;

        va_ngbdisp[20]=dxplus+dyplus+dzminus ;

        va_ngbdisp[21]=dxplus+dyminus+dzplus ;

        va_ngbdisp[22]=dxminus+dyplus+dzplus ;

        va_ngbdisp[23]=dxplus+dyminus+dzminus ;

        va_ngbdisp[24]=dxminus+dyplus+dzminus ;

        va_ngbdisp[25]=dxminus+dyminus+dzplus ;

        va_ngbdisp[26]=dxminus+dyminus+dzminus ;

        va_ngbdisp[0]=0 ;

        for (ngbno=1;ngbno<27;++ngbno) {
            a=cell[j+va_ngbdisp[ngbno]];

            if (a<0) continue;

            if (atom_stat[a]==0) unmoved_ngb[n_unmoved++]=a;
            else if (atom_stat[a]==1) rigid_ngb[n_rigid++]=a;
            else if (atom_stat[a]==2) flxbl_ngb[n_flxbl++]=a;

            while ((a=pnt[a])>=0) {
                if (atom_stat[a]==0) unmoved_ngb[n_unmoved++]=a;
                else if (atom_stat[a]==1) rigid_ngb[n_rigid++]=a;
                else if (atom_stat[a]==2) flxbl_ngb[n_flxbl++]=a;
            }
        }

        //cout<<"select neighbours : found "<<n_unmoved<<" unmoved neighbours\n"
        //<<n_rigid<<" rigid neighbours and "<<n_flxbl<<" flexible neighbours\n"
        //<<"for atom "<<i<<"\n";
        //cout <<"cell id for atom "<<i<<" is "<<cell[atomloc[i]]<<": ("<<lxyz[3*i]<<", "<<lxyz[3*i+1]<<", "<<lxyz[3*i+2]<<")\n";
        //cout <<"and it's neighbours are ...\n";
        //for (a=0;a<n_unmoved;++a) {
        //int ja=unmoved_ngb[a];
        //int jb=atomloc[ja];
        //cout <<"unmoved "<<ja<<" : "<<jb<<": ("<<lxyz[3*ja]<<", "<<lxyz[3*ja+1]<<", "<<lxyz[3*ja+2]<<")\n";
        //}
        //for (a=0;a<n_rigid;++a) {
        //int ja=rigid_ngb[a];
        //int jb=atomloc[ja];
        //cout <<"rigid "<<ja<<" : "<<jb<<": ("<<lxyz[3*ja]<<", "<<lxyz[3*ja+1]<<", "<<lxyz[3*ja+2]<<")\n";
        //}
        //for (a=0;a<n_flxbl;++a) {
        //int ja=flxbl_ngb[a];
        //int jb=atomloc[ja];
        //cout <<"flexible "<<ja<<" : "<<jb<<": ("<<lxyz[3*ja]<<", "<<lxyz[3*ja+1]<<", "<<lxyz[3*ja+2]<<")\n";
        //}
        //getchar();
        if (n_unmoved>=MAXNGB||n_rigid>=MAXNGB||n_flxbl>=MAXNGB) {
            prf::cerr<<"ERROR! Maximum expected number of neighbours ="
            <<MAXNGB<<"\n"
            <<"But n_unmoved = "<<n_unmoved<<", n_flxbl="<<n_flxbl
            <<", n_rigid="<<n_rigid<<"\n";
            exit(1);
        }
    }

    int ExVol::checksubsets(Update *updt)
    {
        int i,j;
        int moved_st, moved_end,rigid0,rigid1;
        double ds2,distold,distnew;
        Vector3 dxv;
        int nerrors=0;

        for (i=0;i<NTO;++i) if (atom_stat[i]!=0) {
                prf::cout<<"broke update start detection loop at "<<i
                <<" with atom status "<<atom_stat[i]<<"\n" ;
                break;
            }

        moved_st=i;

        for (;i<NTO;++i) if (atom_stat[i]==0) {
                prf::cout<<"broke update end detection loop at "<<i
                <<" with atom status "<<atom_stat[i]<<"\n";
                break;
            }

        moved_end=i;

        for (i=0;i<NTO;++i) {
            if (i>=moved_st && i<moved_end) continue;

            dxv=AtomCoordinates::diff_from_backup(i);

            ds2=dxv.dot(dxv);

            if (ds2>1e-7) {
                prf::cout<<"the atom "<<i
                <<" is marked as unmoved, but it's distance "
                <<"from the backup coordinates is "<<ds2<<"\n";
                ++nerrors;
            }
        }

        for (i=moved_st;i<moved_end;++i) {
            dxv=AtomCoordinates::diff_from_backup(i);
            ds2=dxv.dot(dxv);

            if (ds2<1e-7) {
                prf::cout<<"the atom "<<i
                <<" is marked as moved under update "<<updt->Name()
                <<", but it's distance from the backup coordinates is "
                <<ds2<<"\n";
                ++nerrors;
            }
        }

        for (i=moved_st;i<moved_end;++i) {
            if (atom_stat[i]==1) break;
        }

        rigid0=i;

        for (;i<moved_end;++i) {
            if (atom_stat[i]!=1) break;
        }

        rigid1=i;

        prf::cout<<"atoms "<<rigid0<<" through "<<rigid1
        <<" have been marked as the rigid subset\n";

        for (i=rigid0;i<rigid1;++i) {
            for (j=i+1;j<rigid1;++j) {
                AtomCoordinates::switch_to_backups();
                distold=AtomCoordinates::dist2(i,j);
                AtomCoordinates::switch_to_regular();
                distnew=AtomCoordinates::dist2(i,j);

                if (fabs(distnew-distold)>1e-7) {
                    prf::cout<<"atoms "<<i<<" and "<<j
                    <<" are supposed to be in the rigid part but..\n";
                    prf::cout<< "old distance = "<<distold
                    <<" new distance = "<<distnew<<"\n";
                    ++nerrors;
                    //getchar();
                }
            }
        }

        prf::cout<<"checksubsets: number of errors for update "

        <<updt->Name()<<" = "<<nerrors<<"\n";
        prf::cout<<"update start = "<<moved_st
        <<" update end = "<<moved_end<<"\n";
        return 0;
    }

}

