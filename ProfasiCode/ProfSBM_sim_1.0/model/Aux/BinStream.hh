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

#ifndef BINSTREAM_HH
#define BINSTREAM_HH
#include <vector>
#include "fileutils.hh"

class BinStream
{
public:
    BinStream();
    BinStream(const BinStream &);
    ~BinStream();
    BinStream &operator=(const BinStream &b);
    void data(const char *, unsigned long nbyts);
    void set_data_attribs(std::string endianness, std::vector<unsigned int> si);
    BinStream &operator()(unsigned long pos);
    BinStream &operator>>(char &c);
    BinStream &operator>>(unsigned char &c);
    BinStream &operator>>(int &i);
    BinStream &operator>>(unsigned int &u);
    BinStream &operator>>(long &l);
    BinStream &operator>>(unsigned long &u);
    BinStream &operator>>(float &f);
    BinStream &operator>>(double &d);
    BinStream &operator>>(long double &d);
private:
    void tmpcopy(prf_utils::VarType v);
    int dataendianness,myendianness;
    unsigned long loc, bufsize;
    unsigned short origin;
    std::vector<unsigned int> sizes;
    char var[40],*buf;
};

#endif // BINSTREAM_HH
