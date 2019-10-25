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

#include "ConnectionsMatrix.hh"
#include "profasi_io.hh"

using namespace prf;
using std::valarray;

ConnectionsMatrix::ConnectionsMatrix() : nto(0) {}

ConnectionsMatrix::~ConnectionsMatrix() {}

void ConnectionsMatrix::set_connection(int uid1,int uid2)
{
    if (uid1>=nto|| uid2>=nto) {
        prf::cerr<<"can not create connection for ids "<<uid1<<", "<<uid2<<"\n"
        <<"ids must be in the range 0 to "<<nto<<"\n";
        return;
    }

    if (abs(uid1-uid2)>35) {
        Logger(15)<<"Can not create connection for ids "<<uid1<<", "<<uid2<<" as "
        <<"the ids are separated by more than 35. \n";
        return;
    }

    con[72*uid1+36+uid1-uid2]=true;

    con[72*uid2+36+uid2-uid1]=true;
}

void ConnectionsMatrix::NAtoms(int i)
{
    nto=i;
    con.resize(nto*72,false);
}

