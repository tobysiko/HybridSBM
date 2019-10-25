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

#include "AtomDistance.hh"

using namespace prf_utils;
using std::string;
using std::vector;

/**
* \page opt_AtomDistance AtomDistance
\section options Available options
<b>between</b> <i>between category chain_id atom_index1 atom_index2</i><br>
The "between" option must be followed by 4 fields. The last two values give atom indices. The chain_id field tells which chain the atoms are in, if the indices refer to atom indices in a particular chain. The category field can have 3 possible values: backbone/chain/population. This sets the indices to refer to indices within a category of atoms.
\section examples Examples
new_obs AtomDistance adist backbone 0 0 32 <br>
This creates a new task to measure the distance between the first and 33rd atoms along the backbone of chain 0, and calls the measurement "adist" in log files and the averages file. <br><br>
new_obs AtomDistance adistB chain 0 12 65 <br>
This creates a distance measurement between the 13th and 66th atoms of chain 0, whatever those atoms might be, and calls the measurement "adistB". <br><br>
new_obs AtomDistance adistC population 0 23 256 <br>
This creates a distance measurement between the 24th and 256th atom in the population, irrespective of which chain they are on. The chain_id field is, in this case, meaningless.

\sa prf::AtomDistance
*/

namespace prf
{
    AtomDistance::AtomDistance() : a1(0),a2(0) {Name("Dist");}

    AtomDistance::~AtomDistance() {}

    int AtomDistance::init_obs()
    {
        if (Observable::init_obs()==0) return 0;

        vector<string> parts;

        for (size_t i =0;i<usrcmd.size();++i) {
            parts.clear();
            split(usrcmd[i],parts);

            if (parts[0]==string("between")) {
                if (parts.size()<5) {
                    prf::cerr<<"AtomDistance ("<<Name()
                    <<")> Syntax error in command \"between\". Usage:\n"
                    <<"between backbone 0 0 65\n or \n"
                    <<"between chain 5 0 289\n or \n"
                    <<"between population 0 0 654\n\n"
                    <<"The first number in each case means the chain id"
                    <<"which is ignored if the category is population.\n"
                    <<"The other two numbers are atom indices in "
                    <<"that particular category.\n";
                    continue;
                }

                int i1,i2,i3,i4,ich;

                ich=atoi(parts[2].c_str());
                i1=atoi(parts[3].c_str());
                i2=atoi(parts[4].c_str());

                if (i1<0) i1=0;

                if (i2<0) i2=0;

                if (parts[1]==string("backbone")) {
                    if (i1>=(i3=3*p->Chain(ich)->numAminoAcids()))
                        i1=i3-1;

                    if (i2>=(i3=3*p->Chain(ich)->numAminoAcids()))
                        i2=i3-1;

                    i3=p->Chain(ich)->backbone()->atom(i1).UniqueId();

                    i4=p->Chain(ich)->backbone()->atom(i2).UniqueId();
                } else if (parts[1]==string("chain")) {
                    if (i1>=(i3=p->Chain(ich)->numberOfAtoms())) i1=i3-1;

                    if (i2>=(i3=p->Chain(ich)->numberOfAtoms())) i2=i3-1;

                    i3=p->Chain(ich)->begin_atom() +i1;

                    i4=p->Chain(ich)->begin_atom() +i2;
                } else {
                    if (i1>=(i3=p->NumberOfAtoms())) i1=i3-1;

                    if (i2>=(i3=p->NumberOfAtoms())) i2=i3-1;

                    i3=i1;

                    i4=i2;
                }

                set_atoms(i3,i4);

                Logger()(log_thres)<<Name()
                <<"> All set to measure the distance between "
                <<"atoms "<<i3<<" and "<<i4<<" corresponding to "
                <<parts[1]<<" "<<i1<<", "<<i2<<"\n";
            }
        }

        return 1;
    }

    double AtomDistance::evaluate()
    {
        return AtomCoordinates::dist(a1,a2);
    }

    void AtomDistance::rangeEstimate(double &xmin,double &xmax)
    {
        xmin=0;xmax=AtomCoordinates::boxL();
        if (!userbinsz) xbin0=(xmax0-xmin0)/nbins0;
    }

}
