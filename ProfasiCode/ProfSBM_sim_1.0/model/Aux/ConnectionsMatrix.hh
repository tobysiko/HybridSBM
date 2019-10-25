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

#ifndef ConnectionsMatrix_HH
#define ConnectionsMatrix_HH
#include <iostream>
#include <valarray>
#include <bitset>

namespace prf
{

    class ConnectionsMatrix
    {
    public:
        ConnectionsMatrix();
        ~ConnectionsMatrix();
        inline bool operator()(int uid1,int uid2) const
        {if((uid2-=uid1)>35 || uid2<-35) return false; else return con[72*uid1+36-uid2];}

        void NAtoms(int i);
        void set_connection(int uid1,int uid2);
        void write_aa_map(int);
    private:
        std::valarray<bool> con;
        int nto;
    };
}

#endif
