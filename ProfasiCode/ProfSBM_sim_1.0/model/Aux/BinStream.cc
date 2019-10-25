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

#include "BinStream.hh"

using namespace prf_utils;

BinStream::BinStream()
{
    dataendianness=myendianness=prf_utils::running_on_little_endian()?0:1;
    if (myendianness==1) origin=20; else origin=0;
    bufsize=loc=0;
    sizes.clear();
    sizes.resize(7,0);
    sizes[CHAR]=sizeof(char);
    sizes[INT]=sizeof(int);
    sizes[LONG]=sizeof(long);
    sizes[FLOAT]=sizeof(float);
    sizes[DOUBLE]=sizeof(double);
    sizes[LONG_DOUBLE]=sizeof(long double);
    sizes[UNKNOWN]=0;
    buf=0;
}

BinStream::~BinStream()
{
    if (buf) delete(buf);
}

BinStream::BinStream(const BinStream &b)
{
    myendianness=prf_utils::running_on_little_endian()?0:1;
    dataendianness=b.dataendianness;
    if (myendianness==1) origin=20; else origin=0;
    bufsize=b.bufsize;
    loc=b.loc;
    sizes=b.sizes;
    if (bufsize!=0) {
        buf=new char[bufsize];
    }
}

BinStream & BinStream::operator=(const BinStream &b)
{
    if (this!=&b) {
        myendianness=prf_utils::running_on_little_endian()?0:1;
        dataendianness=b.dataendianness;
        if (myendianness==1) origin=20; else origin=0;
        bufsize=b.bufsize;
        loc=b.loc;
        sizes=b.sizes;
        if (bufsize!=0) {
            buf=new char[bufsize];
        }
    }
    return *this;
}

void BinStream::data(const char *d, unsigned long nbyts)
{
    if (nbyts!=bufsize) {
        if (buf) delete buf;
        buf=0;
        if (nbyts!=0) buf=new char[nbyts];
        bufsize=nbyts;
    }
    for (unsigned long i=0;i<nbyts;++i) buf[i]=d[i];
    loc=0;
}

void BinStream::set_data_attribs(std::string endianness,
                                 std::vector<unsigned int> si)
{
    if (endianness=="little_endian") dataendianness=0;
    else if (endianness=="big_endian") dataendianness=1;
    sizes[CHAR]=si[CHAR];
    sizes[INT]=si[INT];
    sizes[LONG]=si[LONG];
    sizes[FLOAT]=si[FLOAT];
    sizes[DOUBLE]=si[DOUBLE];
    sizes[LONG_DOUBLE]=si[LONG_DOUBLE];
}

BinStream &BinStream::operator()(unsigned long pos)
{
    if (pos<bufsize) loc=pos; else loc=bufsize;
    return *this;
}

BinStream &BinStream::operator >>(char &c)
{
    if (loc<bufsize) c=buf[loc++];
    return *this;
}

BinStream &BinStream::operator >>(unsigned char &c)
{
    if (loc<bufsize) c=(unsigned char)buf[loc++];
    return *this;
}

void BinStream::tmpcopy(VarType v)
{
    for (int j=0;j<40;++j) var[j]=0;
    if (dataendianness==myendianness) {
        for (unsigned u=0;u<sizes[v];++u)
            var[origin+u]=buf[loc++];
    } else {
        for (int u=(sizes[v]-1);u>=0;--u)
            var[origin+u]=buf[loc++];
    }
}

BinStream &BinStream::operator >>(int &i)
{
    if ((loc+sizes[INT])<=bufsize) {
        tmpcopy(INT);
        if (myendianness==0) {
            i=*(int *) var;
        } else {
            i=*(int *)(var+origin-sizeof(int)+sizes[INT]);
        }
    }
    return *this;
}

BinStream &BinStream::operator >>(long &i)
{
    if ((loc+sizes[LONG])<=bufsize) {
        tmpcopy(LONG);
        if (myendianness==0) {
            i=*(long *) var;
        } else {
            i=*(long *)(var+origin-sizeof(long)+sizes[LONG]);
        }
    }
    return *this;
}

BinStream &BinStream::operator >>(float &i)
{
    if ((loc+sizes[FLOAT])<=bufsize) {
        tmpcopy(FLOAT);
        if (myendianness==0) {
            i=*(float *) var;
        } else {
            i=*(float *)(var+origin-sizeof(float)+sizes[FLOAT]);
        }
    }
    return *this;
}

BinStream &BinStream::operator >>(double &i)
{
    if ((loc+sizes[DOUBLE])<=bufsize) {
        tmpcopy(DOUBLE);
        if (myendianness==0) {
            i=*(double *) var;
        } else {
            i=*(double *)(var+origin-sizeof(double)+sizes[DOUBLE]);
        }
    }
    return *this;
}

BinStream &BinStream::operator >>(long double &i)
{
    if ((loc+sizes[LONG_DOUBLE])<=bufsize) {
        tmpcopy(LONG_DOUBLE);
        if (myendianness==0) {
            i=*(long double *) var;
        } else {
            i=*(long double *)(var+origin-sizeof(long double)+
                               sizes[LONG_DOUBLE]);
        }
    }
    return *this;
}

BinStream &BinStream::operator >>(unsigned int &i)
{
    int j=0;
    operator>>(j);
    i=(unsigned int) j;
    return *this;
}

BinStream &BinStream::operator >>(unsigned long &i)
{
    long j=0;
    operator>>(j);
    i=(unsigned long) j;
    return *this;
}
