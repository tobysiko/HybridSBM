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

#include <sstream>
#include "His2D.hh"
#include "InstructionString.hh"
using prf::cerr;
using std::string;
using std::swap;
using std::ofstream;
using std::vector;

using namespace prf_utils;

His2D::His2D() : xmin(0),xmax(1),ymin(0),ymax(1),xbin(0.1),ybin(0.1),
        nentries(0),nentries_in(0),nxbin(10),nybin(10) {}

His2D::~His2D() {}

His2D::His2D(const His2D &hs) : xmin(hs.xmin),xmax(hs.xmax),
        ymin(hs.ymin),ymax(hs.ymax),
        xbin(hs.xbin),ybin(hs.ybin),
        nentries(hs.nentries),
        nentries_in(hs.nentries_in),
        nxbin(hs.nxbin),nybin(hs.nybin)
{
    xy.allocate(nxbin,nybin);
    xy=hs.xy;
    X=hs.X;
    Y=hs.Y;
}

His2D & His2D::operator=(const His2D &hs)
{
    if (this!=&hs) {
        xmin=hs.xmin;xmax=hs.xmax;
        ymin=hs.ymin;ymax=hs.ymax;
        xbin=hs.xbin;ybin=hs.ybin;
        nxbin=hs.nxbin;nybin=hs.nybin;
        nentries=hs.nentries;
        nentries_in=hs.nentries_in;
        xy.allocate(nxbin,nybin);
        xy=hs.xy;
        X=hs.X;
        Y=hs.Y;
    }

    return *this;
}

void His2D::init()
{
    if ((nxbin*nybin)==0) {
        prf::cerr<<"Warning: His2D of size "<<nxbin<<" X "<<nybin
        <<" requested. Using "
        <<"size 10 X 10 instead \n";
        nxbin=nybin=10;
    }

    xy.allocate(nxbin,nybin);

    if (xmax<xmin) swap(xmax,xmin);

    if (ymax<ymin) swap(ymax,ymin);

    xbin=(xmax-xmin)/nxbin;ybin=(ymax-ymin)/nybin;

    X.resize(nxbin,0);
    Y.resize(nybin,0);
    for (int i=0;i<nxbin;++i) X[i]=xmin+(0.5+i)*xbin;
    for (int i=0;i<nybin;++i) Y[i]=ymin+(0.5+i)*ybin;
    reset();
}

void His2D::reset()
{
    xy*=0;nentries=0;nentries_in=0;
}

int His2D::put(double x, double y)
{
    ++nentries;
    int ix=(int)floor((x-xmin)/xbin);
    int iy=(int)floor((y-ymin)/ybin);

    if (0<=ix && ix<nxbin && 0<=iy && iy<nybin) {
        xy[ix][iy]+=1;nentries_in++;
        return 1;
    } else return 0;
}

double His2D::normalize()
{
    if (nentries!=0) xy*=(1.0/nentries);
    else prf::cerr <<"can't normalize His2D with zero entries\n";

    return nentries;
}

His2D & His2D::operator+=(const His2D &h1)
{
    int nerr=0;

    if (nxbin!=h1.nxbin) {
        prf::cerr<<"can't add  His2D with different number of xbins\n";
        ++nerr;
    }

    if (nybin!=h1.nybin) {
        prf::cerr<<"can't add His2D with different number of ybins\n";
        ++nerr;
    }

    double xspread=(xmax-xmin),yspread=(ymax-ymin);

    if (fabs(xmin-h1.xmin)>xspread/1000) {
        prf::cerr<<"can't add His2D with different xmin\n";
        ++nerr;
    }

    if (fabs(xmax-h1.xmax)>xspread/1000) {
        prf::cerr<<"can't add His2D with different xmax\n";
        ++nerr;
    }

    if (fabs(ymin-h1.ymin)>yspread/1000) {
        prf::cerr<<"can't add His2D with different ymin\n";
        ++nerr;
    }

    if (fabs(ymax-h1.ymax)>yspread/1000) {
        prf::cerr<<"can't add His2D with different ymax\n";
        ++nerr;
    }

    if (nerr==0) {
        xy+=h1.xy;
        nentries+=h1.nentries;nentries_in+=h1.nentries_in;
    } else {
        prf::cerr<<"addition failed because of "<<nerr<<" errors\n";
    }

    return *this;
}

void His2D::Export(const char *filename, int stl)
{
    switch (stl) {
        case 1: {
            ofstream fout(filename);
            for (int i=0;i<nxbin;++i) {
                for (int j=0;j<nybin;++j) {
                    fout<<xval(i)<<"   "<<yval(j)<<"   "<<xy[i][j]<<"\n";
                    //fout<<xval(i)<<"   "<<yval(j)<<"   "<<xy[j][i]<<"\n";
                }

                fout<<"\n";
            }

            fout.close();
            break;
        }
        case 2: {
            string xfile=string(filename)+"_x";
            string yfile=string(filename)+"_y";
            string xyfile=string(filename);
            ofstream fout(xfile.c_str());

            for (int i=0;i<nxbin;++i) fout<<xval(i)<<"\n";

            fout.close();fout.open(yfile.c_str());

            for (int i=0;i<nybin;++i) fout<<yval(i)<<"\n";

            fout.close();fout.open(xyfile.c_str());

            for (int i=0;i<nxbin;++i) {
                for (int j=0;j<nybin;++j) {
                    fout<<xy[i][j]<<"  ";
                }

                fout<<"\n";
            }

            fout.close();
            break;
        }
    case 3:
    default: {
            std::ofstream fout(filename);
            fout << "# data_layout profasi_h2d_v3\n";
            fout << "# nxbin "<<nxbin<<"\n";
            fout << "# nybin "<<nybin<<"\n";
            fout << "# nentries "<<nentries<<"\n";
            fout << "# nentries_in "<<nentries_in<<"\n";
            fout << "# xmin "<<xmin<<"\n";
            fout << "# xmax "<<xmax<<"\n";
            fout << "# ymin "<<ymin<<"\n";
            fout << "# ymax "<<ymax<<"\n";
            fout <<"0 ";
            for (int j=0;j<nybin;++j) fout << yval(j)<<" ";
            fout <<"\n";
            for (int i=0;i<nxbin;++i) {
                fout<<xval(i)<<" ";
                for (int j=0;j<nybin;++j) {
                    fout<<xy[i][j]<<" ";
                }
                fout<<"\n";
            }
            fout.close();
        }
    };
}

int His2D::Import(const char *filename)
{
    if (STestFile_r(filename)==0) {
        prf::cerr<<"His2D> Unable to open "<<filename<<" for importing data.\n";
        return 0;
    }
    std::ifstream fin(filename);
    std::string line;
    getline(fin,line);
    line=prf_utils::trim_str(line);
    if (line!="# data_layout profasi_h2d_v3") {
        prf::cerr<<"His2D> Import is only possible for His2D data saved in "
                <<"format profasi_h2d_v3.\n";
        return 0;
    }
    do {
        if (line[0]=='#') {
            InstructionString cmds(line.substr(2));
            if (cmds.head()=="nxbin") nxbin=atoi(cmds.tail().str().c_str());
            else if (cmds.head()=="nybin") nybin=atoi(cmds.tail().str().c_str());
            else if (cmds.head()=="nentries") nentries=atoi(cmds.tail().str().c_str());
            else if (cmds.head()=="nentries_in") nentries_in=atoi(cmds.tail().str().c_str());
            else if (cmds.head()=="xmin") xmin=strtod(cmds.tail().str().c_str(),NULL);
            else if (cmds.head()=="xmax") xmax=strtod(cmds.tail().str().c_str(),NULL);
            else if (cmds.head()=="ymin") ymin=strtod(cmds.tail().str().c_str(),NULL);
            else if (cmds.head()=="ymax") ymax=strtod(cmds.tail().str().c_str(),NULL);
        } else break;
    } while (getline(fin,line));

    long nentries_bkp=nentries,nentries_in_bkp=nentries_in;
    init();


    std::istringstream ssin(line);

    double dummy=0;
    ssin>>dummy;

    for (int j=0;j<nybin;++j) ssin>>Y[j];

    for (int i=0;i<nxbin;++i) {
        fin>>X[i];
        for (int j=0;j<nybin;++j) {
            fin>>xy[i][j];
        }
    }
    nentries=nentries_bkp;
    nentries_in=nentries_in_bkp;
    fin.close();
    return 1;
}

void His2D::Projection_X(vector<double> &hsx)
{
    hsx.resize(NXbins(),0.0);

    for (int i=0;i<NXbins();++i) {
        double ysum=0;

        for (int j=0;j<NYbins();++j) {
            ysum+=xy[i][j];
        }

        hsx[i]=ysum;
    }
}

void His2D::Projection_Y(vector<double> &hsy)
{
    hsy.resize(NYbins(),0.0);

    for (int i=0;i<NYbins();++i) {
        double xsum=0;

        for (int j=0;j<NXbins();++j) {
            xsum+=xy[j][i];
        }

        hsy[i]=xsum;
    }
}

void His2D::save_state(const char *filename)
{
    FILE *fp;
    fp=fopen(filename,"w");
    fwrite(&nxbin,sizeof(int),1,fp);
    fwrite(&nybin,sizeof(int),1,fp);
    fwrite(&nentries,sizeof(long),1,fp);
    fwrite(&nentries_in,sizeof(long),1,fp);
    fwrite(&xmin,sizeof(double),1,fp);
    fwrite(&xmax,sizeof(double),1,fp);
    fwrite(&ymin,sizeof(double),1,fp);
    fwrite(&ymax,sizeof(double),1,fp);
    xy.Write(fp);
    fclose(fp);
}

void His2D::read_state(const char *filename)
{
    long n1,n2;
    FILE *fp;
    fp=fopen(filename,"r");
    fread(&nxbin,sizeof(int),1,fp);
    fread(&nybin,sizeof(int),1,fp);
    fread(&n1,sizeof(long),1,fp);
    fread(&n2,sizeof(long),1,fp);
    fread(&xmin,sizeof(double),1,fp);
    fread(&xmax,sizeof(double),1,fp);
    fread(&ymin,sizeof(double),1,fp);
    fread(&ymax,sizeof(double),1,fp);
    init();
    nentries=n1;nentries_in=n2;
    xy.Read(fp);
    fclose(fp);
}
