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

#include "GroupProps.hh"
#include "../Aux/profasi_io.hh"
#include "../Aux/fileutils.hh"
#include <algorithm>

using std::string;

namespace prf
{
    GroupProps::GroupProps()
    {
        isaa=iseg=issynth=isnat=hyph=chgd=false;
        natoms=nhvatoms=nsdof=0;
        comname="VoidGroup";threelet="ERR";
    }

    GroupProps::GroupProps(string lbls)
    {
        isaa=iseg=issynth=isnat=false;
        natoms=nhvatoms=nsdof=0;
        comname="VoidGroup";threelet="ERR";
        init(lbls);
    }

    GroupProps::~GroupProps() {}

    void GroupProps::init(string lbls)
    {
        Logger blog;
        blog(150)<<comname<<":\n";
        thelabel.clear();
        natoms=prf_utils::split_str<std::vector<string> >(lbls,',',thelabel);
        nhvatoms=0;

        for (int i=0;i<natoms;++i) {
            string nm=thelabel[i];
            blog(150)<<"atom "<<i<<" = "<<nm<<"\n";
            theindex[nm]=i+1;

            if (nm[1]!='H') ++nhvatoms;
        }

        nsdof=0;
    }

    void GroupProps::torsion_dof(int ai[4])
    {
        Logger blog;
        defatoms.insert(defatoms.end(),ai,ai+4);
        ++nsdof;
        blog(150)<<"Incremented side chain degrees of freedom to "<<nsdof<<"\n";
    }

    void GroupProps::torsion_dof(int a1, int a2, int a3, int a4)
    {
        Logger blog;
        defatoms.push_back(a1);defatoms.push_back(a2);
        defatoms.push_back(a3);defatoms.push_back(a4);
        ++nsdof;
        blog(150)<<"Incremented side chain degrees of freedom to "<<nsdof<<"\n";
    }

    void GroupProps::dof_def_atoms(int idof, int ai[4])
    {
        for (int i=0;i<4;++i) ai[i]=defatoms[4*idof+i];
    }

    void GroupProps::dof_def_atoms(int idof, int &a1,int &a2, int &a3,int &a4)
    {
        a1=defatoms[4*idof]; a2=defatoms[4*idof+1];
        a3=defatoms[4*idof+2]; a4=defatoms[4*idof+3];
    }

    void GroupProps::set_type(string typ)
    {
        if (typ.find("aminoacid")<typ.size()) isaa=true;

        if (typ.find("natural")<typ.size()) isnat=true;

        if (typ.find("synthetic")<typ.size()) issynth=true;

        if (typ.find("endgroup")<typ.size()) iseg=true;

        if (typ.find("hydrophobic")<typ.size()) hyph=true;

        if (typ.find("charged")<typ.size()) chgd=true;

        typedescr=typ;

        if (issynth) isnat=false;
    }

    AtomKind GroupProps::species(int i)
    {
        switch (thelabel[i][1]) {
            case 'C': return carbon;
            case 'N': return nitrogen;
            case 'O': return oxygen;
            case 'S': return sulfur;
            case 'H':
            default: return hydrogen;
        };
    }

    bool GroupProps::valid_label(std::string lbl) const
    {
        return (std::find(thelabel.begin(),thelabel.end(),lbl)!=thelabel.end());
    }
}
