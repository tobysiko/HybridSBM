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

#include "Permutation.hh"

using std::vector;

namespace prf_utils
{
    Permutation::Permutation() : nmembers(0) {}

    Permutation::Permutation(int n) : nmembers(n),valperm(n,0),trperm(n,0) {}

    Permutation::~Permutation() {}

    Permutation & Permutation::operator=(Permutation & prm)
    {
        if (this!=&prm) {
            nmembers=prm.nmembers;
            valperm=prm.valperm;
            trperm=prm.trperm;
        }

        return (*this);
    }

    bool Permutation::IsIdentity()
    {
        bool ans=true;

        for (int i=0;i<nmembers;++i) {
            if (valperm[i]!=i) {ans=false;break;}
        }

        return ans;
    }

    void Permutation::reset(int siz)
    {
        nmembers=siz;
        valperm=vector<int>(nmembers);
        trperm=vector<int>(nmembers);
    }

    void Permutation::SetToIdentity()
    {
        for (int i=0;i<nmembers;++i) {
            valperm[i]=i;
            trperm[i]=i;
        }
    }

    bool Permutation::Element(int i)
    {
        if (0<=i && i<nmembers) return true; else return false;
    }

    int Permutation::Assign(int indx, int val)
    {
        if (Element(indx)) {
            valperm[indx]=val;
            trperm[val]=indx;
            return 1;
        } else return 0;
    }

    int Permutation::FlipContents(int i, int j)
    {
        if (Element(i)&&Element(j)) {
            int tmpint=valperm[i];
            valperm[i]=valperm[j];
            valperm[j]=tmpint;
            trperm[valperm[j]]=j;
            trperm[valperm[i]]=i;
            return 1;
        } else return 0;
    }

    int Permutation::FlipLocations(int i, int j)
    {
        int tmpint=trperm[j];
        trperm[j]=trperm[i];
        trperm[i]=tmpint;
        valperm[trperm[i]]=i;
        valperm[trperm[j]]=j;
        return 1;
    }

    prf::Output & operator<<(prf::Output & os,
                             prf_utils::Permutation &prm)
    {
        for (int i=0;i<prm.nmembers;++i) {
            os<<prm[i]<<"\t";

            if ((1+i)%10==0) os <<"\n";
        }

        os<<"\n";

        return os;
    }
}
