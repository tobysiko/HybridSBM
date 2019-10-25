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

#include "His1D.hh"
#include "profasi_io.hh"
#include "fileutils.hh"
#include <cstdio>
#include <sstream>
#include <deque>
#include <fstream>

using prf::cerr;
using std::swap;
using std::vector;

using namespace prf_utils;
His1D::His1D() : xmin(0),xmax(1),xbin(0.1),nbin(100),nblk(1),
    nam("unnamed") {}

His1D::His1D(int numblocks) : xmin(0),xmax(1),xbin(0.1),nbin(100),
        nblk(numblocks),nam("unnamed") {}

His1D::His1D(double xmn,double xmx,int npnts,int numblocks) :
        xmin(xmn),xmax(xmx),xbin(0.1),nbin(npnts),nblk(numblocks),
        nam("unnamed") {}

His1D::His1D(double xmn, double xmx, double bnsz, int numblocks) :
        xmin(xmn),xmax(xmx),xbin(bnsz),nblk(numblocks),nam("unnamed")
{
    set_bin_size(bnsz);
}

His1D::His1D(const His1D &hs) :

        xmin(hs.xmin),xmax(hs.xmax),xbin(hs.xbin),nbin(hs.nbin),nblk(hs.nblk)
{
    nentries=hs.nentries;
    nentries_in=hs.nentries_in;
    xy.allocate(nblk,nbin);
    xy=hs.xy;
    nam=hs.nam;
}

His1D & His1D::operator=(const His1D &hs)
{
    if (this!=&hs) {
        xmin=hs.xmin;
        xmax=hs.xmax;
        xbin=hs.xbin;
        nbin=hs.nbin;
        nblk=hs.nblk;
        nentries=hs.nentries;
        nentries_in=hs.nentries_in;
        xy.allocate(nblk,nbin);
        xy=hs.xy;
        nam=hs.nam;
    }

    return *this;
}

His1D::~His1D() {}

void His1D::set_bin_size(double sz)
{
    xbin=sz;
    if (xmax<xmin) swap(xmin,xmax);
    nbin=ceil((xmax-xmin)/xbin);
    double extn=nbin*xbin-(xmax-xmin);
    xmax+=0.5*extn;
    xmin-=0.5*extn;
}

void His1D::init()
{
    if (nbin<=0) {
        prf::cerr<<"Warning: His1D of size "<<nbin
        <<" requested. Using size 100 instead \n";
        nbin=100;
    }

    xy.allocate(nblk,nbin);

    if (xmax<xmin) swap(xmax,xmin);

    xbin=(xmax-xmin)/nbin;

    reset();
}

void His1D::reset()
{
    xy*=0;nentries_in=vector<double>(nblk,0);nentries=vector<double>(nblk,0);
}

int His1D::put(double x, int iblk)
{
    nentries[iblk]+=1;
    //Bug fix in v. 1.0.2. Thanks to : Simon Mitternacht
    //The bug caused a segmentation fault in extremely rare cases where
    //x is rather close to xmax. The condition x<xmax then in some machines
    //did not always ensure that (x-xmin)/xbin < nbin. In perfect machines
    //with infinite precision this should never happen, but since it does,
    //we use this rearrangement.
    int idx= (int) floor((x-xmin)/xbin);

    if (0<=idx&&idx<nbin) {
        xy[iblk][idx]+=1;
        nentries_in[iblk]+=1;
        return 1;
    } else return 0;
}

int His1D::nput(double howmanytimes, double x, int iblk)
{
    nentries[iblk]+=howmanytimes;
    int idx= (int) floor((x-xmin)/xbin);

    if (0<=idx&&idx<nbin) {
        xy[iblk][idx]+=howmanytimes;
        nentries_in[iblk]+=howmanytimes;
        return 1;
    } else return 0;
}

double His1D::unnormalize()
{
    double totnorm=0;

    for (int i=0;i<nblk;++i) {
        totnorm=0;

        for (int j=0;j<nbin;++j) {
            totnorm+=xy[i][j];
        }

        if (fabs(totnorm)>1e-7 && fabs(totnorm-nentries_in[i])>1e-7) {
            for (int j=0;j<nbin;++j) {
                xy[i][j]*=(nentries[i]);
            }
        }
    }

    return totnorm;
}

double His1D::normalize()
{
    double totnorm=0,blknorm=1;

    for (int i=0;i<nblk;++i) {
        if (nentries[i]!=0) {
            blknorm=1.0/nentries[i];

            for (int j=0;j<nbin;++j) {
                xy[i][j]*=blknorm;
            }
        } else prf::cerr <<"His1D: can't normalize block "<<i
            <<" which has zero entries\n";

        totnorm+=nentries[i];
    }

    return totnorm;
}

His1D & His1D::operator+=(His1D &h1)
{
    int nerr=0;

    if (nbin!=h1.nbin) {
        prf::cerr<<"can't add  His1D with different number of bins\n";
        ++nerr;
    }

    if (nblk!=h1.nblk) {
        prf::cerr<<"can't add His1D with different number of blocks\n";
        ++nerr;
    }

    double xspread=(xmax-xmin);

    if (fabs(xmin-h1.xmin)>xspread/1000) {
        prf::cerr<<"can't add His1D with different xmin\n";
        ++nerr;
    }

    if (fabs(xmax-h1.xmax)>xspread/1000) {
        prf::cerr<<"can't add His1D with different xmax\n";
        ++nerr;
    }

    if (nerr==0) {
        xy+=h1.xy;

        for (int i=0;i<nblk;++i) {
            nentries[i]+=h1.nentries[i];nentries_in[i]+=h1.nentries_in[i];
        }
    } else {
        prf::cerr<<"addition failed because of "<<nerr<<" errors\n";
    }

    return *this;
}

void His1D::Export(const char *filename, int normmode, int datlayout)
{
    FILE *fp=fopen(filename,"w");
    fprintf(fp,"# data_layout profasi_his_v%d\n",datlayout);
    fprintf(fp,"# nblocks  %d\n",nblk);
    fprintf(fp,"# nbins  %d\n",nbin);
    fprintf(fp,"# xmin  %.16f\n",xmin);
    fprintf(fp,"# xmax  %.16f\n",xmax);
    fprintf(fp,"# normalization_mode %d\n",normmode);

    for (int i=0;i<nblk;++i) {
        fprintf(fp,"# block_entries %d %f %f\n",i,nentries[i],
                nentries_in[i]);
    }
    std::ostringstream ost;
    ost.precision(16);
    ost.setf(std::ios_base::scientific);

    if (datlayout==1 && nblk!=1) {
        for (int iblk=0;iblk<nblk;++iblk) {
            double blknorm=1.0;
            if (nentries[iblk]!=0) {
                if (normmode==1) blknorm/=nentries[iblk];
                else if (normmode==2) blknorm/=(nentries[iblk]*xbin);
            }

            for (int i=0;i<nbin;++i) {
                ost<<iblk<<"\t"<<xval(i)<<"\t"<<(xy[iblk][i]*blknorm)<<"\n";
            }
            ost<<"\n";
        }
    } else {
        std::vector<double> blknorms(nblk,1.0);
        for (int iblk=0;iblk<nblk;++iblk) {
            if (nentries[iblk]!=0) {
                if (normmode==1) blknorms[iblk]/=nentries[iblk];
                else if (normmode==2) blknorms[iblk]/=(nentries[iblk]*xbin);
            }
        }

        for (int i=0;i<nbin;++i) {
            ost<<xval(i);
            for (int iblk=0;iblk<nblk;++iblk) {
                ost<<"\t"<<(xy[iblk][i]*blknorms[iblk]);
            }
            ost<<"\n";
        }

    }

    fprintf(fp,"%s",ost.str().c_str());
    fclose(fp);
}

int His1D::Import(const char *filename)
{
    double n1,n2;
    int iblk=0,nrmmd=1, datlayout=1;
    std::deque<std::string> lines;
    std::deque<std::deque<std::string> > later;
    std::string line;

    if (TestFile_r(filename)==0) return 0;

    std::ifstream fin(filename);
    while (getline(fin,line)) lines.push_back(trim_str(line));
    fin.close();

    size_t iline=0;
    while (iline<lines.size() && lines[iline][0]=='#') {
        std::string metainfo=lines[iline].substr(1);
        std::deque<std::string> tokens;
        split(metainfo,tokens);
        if (tokens[0]=="data_layout") {
            if (tokens[1]=="profasi_his_v2") datlayout=2;
        } else if (tokens[0]=="nblocks") nblk=atoi(tokens[1].c_str());
        else if (tokens[0]=="nbins") nbin=atoi(tokens[1].c_str());
        else if (tokens[0]=="xmin") xmin=strtod(tokens[1].c_str(),NULL);
        else if (tokens[0]=="xmax") xmax=strtod(tokens[1].c_str(),NULL);
        else if (tokens[0]=="normalization_mode") nrmmd=atoi(tokens[1].c_str());
        else if (tokens[0]=="block_entries") later.push_back(tokens);
        ++iline;
    }

    init();

    for (size_t i=0;i<later.size();++i) {
        iblk=atoi(later[i][1].c_str());
        n1=strtod(later[i][2].c_str(),NULL);
        n2=strtod(later[i][3].c_str(),NULL);
        nentries[iblk]=n1;
        nentries_in[iblk]=n2;
    }

    double dummyx,dummyy;
    if (nblk==1) datlayout=2;

    if (datlayout==2) {
        int ibin=0;
        for (;iline<lines.size();++iline) {
            line=trim_str(lines[iline]);
            if (line.empty()) continue;
            std::istringstream ssin(line);
            ssin>>dummyx;
            for (int j=0;j<nblk;++j) {
                ssin>>dummyy;
                if (nentries[j]!=0) xy[j][ibin]=dummyy;
                else xy[j][ibin]=0;
            }
            ++ibin;
        }
    } else {
        int ibin=0;
        std::deque<std::string> tokens;
        for (;iline<lines.size();++iline) {
            line=trim_str(lines[iline]);
            if (line.empty()) continue;
            tokens.clear();
            split(line,tokens);
            iblk=atoi(tokens[0].c_str());
            if (iblk<0 || iblk>=nblk || tokens.size()<3) continue;
            dummyx=strtod(tokens[1].c_str(),NULL);
            dummyy=strtod(tokens[2].c_str(),NULL);
            ibin=(int) ((dummyx-xmin)/xbin);
            if (ibin<0 || ibin>=nbin) continue;
            if (nentries[iblk]!=0) {
                xy[iblk][ibin]=dummyy;
            } else xy[iblk][ibin]=0;
        }
    }

    std::vector<double> sums(nblk,0);
    for (size_t i=0;i<sums.size();++i) {
        for (int j=0;j<nbin;++j) {
            sums[i]+=xy[i][j];
        }
        if (nentries_in[i]!=0 and sums[i]!=0 and
            fabs(sums[i]-nentries_in[i])/nentries_in[i] > 1e-4) {
            for (int j=0;j<nbin;++j) xy[i][j]*=(nentries_in[i]/sums[i]);
        }
    }

    prf::Logger(10)<<"Import: data layout "<<
            (datlayout==2?"profasi_his_v2":"profasi_his_v1")<<"\n";

    return 1;
}

void His1D::import_data(Matrix<double> &mtx)
{
    if (mtx.dim1()==(size_t) nblk && mtx.dim2()== (size_t) nbin) xy=mtx;
}
