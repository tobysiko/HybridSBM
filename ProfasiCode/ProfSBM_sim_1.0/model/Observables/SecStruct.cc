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

#include "SecStruct.hh"
#include "../Aux/Constants.hh"
#include <fstream>
#include <vector>

using std::ifstream;
using std::vector;
using std::string;

using namespace UnivConstants;

using namespace prf;

using namespace prf_utils;

RCBin::RCBin()
{
    phimin=-pi;phimax=pi;
    psimin=-pi;psimax=pi;
    refvec.clear();
    fixed_his=true;
}

RCBin::~RCBin() {}

void RCBin::his_range(double xmn, double xmx) {}

void RCBin::his_nbins(int n) {}

void RCBin::rangeEstimate(double &x0, double &x1)
{
    if (!refvec.empty()) {
        x0=-0.5/refvec.size();
        xbin0=1.0/refvec.size();
    } else {
        x0=0;
        xbin0=0.01;
    }
    nbins0=refvec.size()+1;
    usernbins=true;
    x1=1-x0;
}

int RCBin::set_res_list(string filename)
{
    if (TestFile_r(filename.c_str())==0) return 0;

    ifstream fin(filename.c_str());

    refvec.clear();

    int dummy;

    while (fin>>dummy)
        if (dummy>=0 && dummy<p->NumberOfResidues())
            refvec.push_back(dummy);
        else prf::cerr<<Name()<<"> Ignoring non-existent residue index "
            <<dummy<<"\n";

    fin.close();

    return 1;
}

double RCBin::evaluate()
{
    double phi=0,psi=0;
    int nin=0;

    for (size_t i=0;i<refvec.size();++i) {
        p->amino_acid(refvec[i])->calcPhiPsi(phi,psi);

        if ((phi>phimin)&&(phi<phimax)
            &&(psi>psimin)&&(psi<psimax)) {
            sttres[i]=1;
            ++nin;
        } else sttres[i]=0;
    }

    return ((double) nin)/refvec.size();
}

/**
* \page opt_RCBin RCBin: Secondary structure analysis based on backbone angles
\section options Available options
<ul>
<li><b>phi_limits</b> phi_limits -1.5708 -3.1416<br>Monitored range in Ramachandran phi angle. </li>
<li><b>psi_limits</b> psi_limits -1.5708 -3.1416<br>Monitored range in Ramachandran psi angle. </li>
<li><b>default_limits</b> default_limits helix<br> Sets up the phi and psi limits as per the default definition of a helix in PROFASI: -90 pi/180 &lt;phi&lt;-30 pi/180 and -77 pi/180 &lt; psi &lt; -17 pi/180. <br> "default_limits strand" will set of the limits according to the default definition of a strand conformation for a residue: -150 pi/180 &lt;phi&lt;-90 pi/180 and 90 pi/180 &lt;psi&lt;150 pi/180. </li>
<li><b>residues_of_chain</b> "residues_of_chain 2" will select all residues of the third chain, barring the two terminal ones for backbone angle monitoring. </li>
<li><b>residues_from_file</b> "residues_from_file index_list_file" will read a list of residues from the given file, and restrict calculations of backbone angle properties to those.</li>
</ul>
\sa prf::RCBin
*/
int RCBin::init_obs()
{
    Logger blog(log_thres);

    if (Observable::init_obs()==0) return 0;

    vector<string> parts;

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);

        if (parts[0]==string("default_limits") && parts.size()>=2) {
            blog<<Name()<<"> Using default limits :";

            if (parts[1]==string("helix") or parts[1]==string("alpha")) {
                phimin=helix_phimin;
                phimax=helix_phimax;
                psimin=helix_psimin;
                psimax=helix_psimax;
                blog<<" (helix) :";
            } else if (parts[1]==string("strand") or parts[1]==string("beta")) {
                phimin=sheet_phimin;
                phimax=sheet_phimax;
                psimin=sheet_psimin;
                psimax=sheet_psimax;
                blog<<" (strand) : ";
            }

            blog<<phimin*radian_in_degrees<<" degrees <= phi < "

            <<phimax*radian_in_degrees<<" degrees,  and "
            <<psimin*radian_in_degrees<<" degrees  <= psi < "
            <<psimax*radian_in_degrees<<" degrees \n";
        }

        if (parts[0]==string("phi_limits") && parts.size()>=3) {
            phimin=strtod(parts[1].c_str(),NULL);
            phimax=strtod(parts[2].c_str(),NULL);
            blog<<Name()<<"> User defined phi range "
            <<phimin<<" to "<<phimax<<"\n";

        }

        if (parts[0]==string("psi_limits") && parts.size()>=3) {
            psimin=strtod(parts[1].c_str(),NULL);
            psimax=strtod(parts[2].c_str(),NULL);
            blog<<Name()<<"> User defined psi range "
            <<psimin<<" to "<<psimax<<"\n";
        }

        if (parts[0]==string("residues_from_file") && parts.size()>=2) {
            set_res_list(parts[1]);
        }

        if (parts[0]==string("residues_of_chain") && parts.size()>=2) {
            int ich=atoi(parts[1].c_str());

            if (ich<0 || ich>=p->NumberOfChains()) {
                prf::cerr<<Name()<<"> Can not take residues from "
                <<"non-existant chain "<<ich<<"\n";
            } else {
                refvec.clear();
                int iaa0=0;

                for (int j=0;j<ich;++j) iaa0+=p->Chain(j)->numAminoAcids();

                for (int j=1;j<p->Chain(ich)->numAminoAcids()-1;++j)
                    refvec.push_back(iaa0+j);

                blog<<Name()<<"> Chose residues 1 to "
                <<p->Chain(ich)->numAminoAcids()-1<<" of chain "
                <<ich<<" for monitoring\n";
            }
        }
    }

    if (refvec.empty()) {
        int iaa0=0;

        for (int ich=0;ich<p->NumberOfChains();++ich) {
            for (int j=1;j<p->Chain(ich)->numAminoAcids()-1;++j) {
                refvec.push_back(iaa0+j);
            }

            iaa0+=p->Chain(ich)->numAminoAcids();
        }
    }

    sttres.resize(refvec.size(),0);

    dndt.resize(ntmp,0);
    prof.allocate(ntmp,refvec.size());
    return 1;
}

void RCBin::avg_fill(int i)
{
    Observable::avg_fill(i);
    if (gathstat) {
        for (size_t j=0;j<refvec.size();++j) prof[i][j]+=sttres[j];
        dndt[i]++;
    }
}

void RCBin::avg_reset()
{
    Observable::avg_reset();
    prof*=0;
    dndt*=0;
}

void RCBin::avg_write(Output &op)
{
    Observable::avg_write(op);
    Output oq((oprefx+Name()+string(".profile")).c_str(),"w");

    for (int i=0;i<ntmp;++i) {
        for (size_t j=0;j<refvec.size();++j) {
            if (dndt[i]>0) oq<<prof[i][j]/dndt[i]<<"  "; else oq<<0.0<<"  ";
        }

        oq<<"\n";
    }

    oq.close();
}
