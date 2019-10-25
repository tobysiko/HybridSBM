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

#include "RCAngleRMSD.hh"
#include "../Aux/Constants.hh"
#include "../Aux/PDBReader.hh"
#include "../Aux/RMSD_Utils.hh"
#include <fstream>

using namespace UnivConstants;

using namespace prf;

using namespace prf_utils;
using std::string;
using std::valarray;
using std::vector;
using std::list;

RCAngleRMSD::RCAngleRMSD()
{
    Name("RCAngleRMSD");
    refres.clear();
    refphi.clear();
    refpsi.clear();
    diff2res.clear();
    p=NULL;
}

RCAngleRMSD::~RCAngleRMSD() {}

void RCAngleRMSD::rangeEstimate(double &x0, double & x1)
{
    x0=0;x1=pi*pi;
    if (!userbinsz) xbin0=0.1;
}

/**
\page opt_RCAngleRMSD RMSD over Ramachandran angles
\section options Available Options
\li \b ref_angle_file ref_angle_file phipsi.list<br> The file with the angles should contain 3 columns: residue index i, phi_i, psi_i. Note that the residue indices not found in the file are not considered for the calculation.
\li \b ref_angle_pdbfile ref_angle_pdbfile molecule.pdb<br> This class can infer the comparison values for the Ramachandran angles from a given structure in the form of a PDB file. Selection rules described in \ref prf_sel_fils apply.

\sa prf::RCAngleRMSD
*/
int RCAngleRMSD::init_obs()
{
    if (Observable::init_obs()==0) return 0;

    vector<string> parts;

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);

        if (parts[0]==string("ref_angle_file") && parts.size()>=2) {
            reference_angle_file(parts[1]);
        }

        if (parts[0]==string("ref_angle_pdbfile") && parts.size()>=2) {
            reference_angle_pdbfile(parts[1]);
        }
    }

    return init();
}

void RCAngleRMSD::reference_angle_list(vector<int> &resind,
                                       vector<double> &phvec,vector<double> &psvec)
{
    Logger blog;
    refres.clear();
    refphi.clear();
    refpsi.clear();
    blog(log_thres)<<"RCAngleRMSD("<<Name()<<")> Listing reference angles\n"
    <<"Residue index \t Phi \t\t Psi\n";

    for (size_t i=0;i<resind.size();++i) {
        if (i>=phvec.size() || i>=psvec.size()) break;

        if (resind[i]<0 || resind[i]>=p->NumberOfResidues()) {
            prf::cerr<<"RCAngleRMSD> residue index "<<resind[i]
            <<" out of range. Ignored.\n";
            continue;
        }

        if (p->amino_acid(resind[i])->hasNTerminal() ||
            p->amino_acid(resind[i])->hasCTerminal()) {
            blog<<"RCAngleRMSD> residue with index "<<resind[i]
            <<"("<<p->amino_acid(resind[i])->Name()
            <<") is at the end of a chain. Ignored.\n";
            continue;
        }

        refres.push_back(resind[i]);

        refphi.push_back(phvec[i]);
        refpsi.push_back(psvec[i]);
        blog<<resind[i]<<"   "<<phvec[i]*radian_in_degrees
        <<"   "<<psvec[i]*radian_in_degrees<<"\n";
    }

    blog<<"RCAngleRMSD> End of listing\n";
}

void RCAngleRMSD::reference_angle_file(string filename)
{
    if (TestFile_r(filename)==0) return;

    std::ifstream fin(filename.c_str());

    int i;

    double ph,ps;

    vector<int> resind;

    vector<double> phvec,psvec;

    while (fin>>i) {
        fin>>ph;fin>>ps;
        resind.push_back(i);
        phvec.push_back(ph);
        psvec.push_back(ps);
    }

    fin.close();

    return reference_angle_list(resind,phvec,psvec);
}

int RCAngleRMSD::reference_angle_pdbfile(string infile)
{
    string fil1;
    string sel1,sel2;
    size_t icolon=infile.find(':');

    if (icolon<infile.size()-1) sel1=string(infile,icolon+1);
    else sel1="1:*";

    sel2="1:*";

    fil1=string(infile,0,icolon);

    PDBReader pdb;

    pdb.set_file(fil1);

    if (pdb.read_matrix()==0) {
        prf::cerr<<"Could not read in information from "<<fil1<<"\n"
        <<"RCAngleRMSD> Inferring comparison angles from PDB file failed.\n";
        return 0;
    }

    list<SelRes> slist1,slist2;

    pdb.mk_selection(sel1,slist1);
    p->mk_selection(sel2,slist2);
    //Sequence alignment, residue listing.
    RMSD_Utils utils;

    if (!utils.seq_match(slist1,slist2)) {
        prf::cerr<<"NativenessQ: Sequences of the two specified objects "
        <<"do not match.\n Performing sequence alignment...\n";
        utils.seq_align(slist1,slist2);

        if (slist1.empty()||slist2.empty()) {
            prf::cerr<<"Attempt to align sequences produced empty lists\n"
            <<"RCAngleRMSD> Inferring comparison angles from PDB file failed.\n";
            return 0;
        }
    }

    list<AtomDescriptor> des1;

    pdb.descriptors(slist1,des1);
    utils.make_filter(string("+BB"));
    utils.apply_filter(des1);
    Shape natv;
    vector<int> vct(des1.size(),0);
    int js2=0;

    for (list<AtomDescriptor>::iterator i=des1.begin();i!=des1.end();++i)  {
        vct[js2++]=i->int_label;
    }

    pdb.export_shape(vct,natv);

    if (natv.NPoints()!=((int) slist2.size() *3)) {
        prf::cerr<<"Exported backbone structure from PDB has "<<natv.NPoints()<<"\n"
        <<"The population wants "<<3*slist2.size()<<" points. \n"
        <<"PDB file "<<infile<<" can not be used for initializing reference "
        <<"angles for RCAngleRMSD\n";
        return 0;
    }

    double ph,ps;

    vector<int> resind;
    vector<double> phvec,psvec;
    Vector3 v0(0,0,0),v1(0,0,0),v2(0,0,0),v3(0,0,0);

    list<SelRes>::iterator it=slist2.begin();

    for (size_t i=0;i<slist2.size();++i,++it) {
        resind.push_back(it->indx_nat);
        v1=natv.Point(3*i+1)-natv.Point(3*i);
        v2=natv.Point(3*i+2)-natv.Point(3*i+1);
        v1.mag(1);v2.mag(1);
        if (i==0) ph=-pi; else ph=v1.torsion(v0,v2);
        if (i==(slist2.size()-1)) ps=pi; else {

            v3=natv.Point(3*i+3)-natv.Point(3*i+2);
            ps=v2.torsion(v1,v3);
            v0=v3;
        }

        phvec.push_back(ph);

        psvec.push_back(ps);
    }

    reference_angle_list(resind,phvec,psvec);

    return 1;
}

int RCAngleRMSD::init()
{
    if (refres.empty()) {
        prf::cerr<<"RCAngleRMSD> Empty reference structure. Initialization failed.\n";
        return 0;
    }

    diff2res.resize(refres.size(),0.0);

    resprops.allocate(ntmp,refres.size());
    dndt.resize(ntmp);
    resprops*=0;
    dndt*=0;
    return 1;
}

double RCAngleRMSD::evaluate()
{
    diff2mean=0;
//Ignore the first and last residues for which the Ramachandran
//angles are not both defined. If there are EndGroups attached,
//these angles would be defined. But for uniformity, we ignore
//them even for that case.

    for (size_t i=0;i<refres.size(); ++i) {
        double tmpdiff1,tmpdiff2;
        AminoAcid *ac=p->amino_acid(refres[i]);
        ac->calcPhiPsi(tmpdiff1,tmpdiff2);
        tmpdiff1-=refphi[i];

        while (tmpdiff1>pi) tmpdiff1-=twoPi;

        while (tmpdiff1<-pi) tmpdiff1+=twoPi;

        tmpdiff2-=refpsi[i];

        while (tmpdiff2>pi) tmpdiff2-=twoPi;

        while (tmpdiff2<-pi) tmpdiff2+=twoPi;

        diff2mean+=tmpdiff1*tmpdiff1+tmpdiff2*tmpdiff2;

        diff2res[i]=tmpdiff1*tmpdiff1+tmpdiff2*tmpdiff2;
    }

    diff2mean/=(2.0* refres.size());

    return sqrt(diff2mean);
}

void RCAngleRMSD::avg_fill(int i)
{
    Observable::avg_fill(i);
    if (gathstat) {
        for (size_t j=0;j<refres.size();++j) resprops[i][j]+=diff2res[j];
        dndt[i]++;
    }
}

void RCAngleRMSD::avg_reset()
{
    Observable::avg_reset();
    resprops*=0;
    dndt*=0;
}

void RCAngleRMSD::avg_write(Output &op)
{
    Observable::avg_write(op);
    Output oq((oprefx+Name()+string(".profile")).c_str(),"w");

    for (int i=0;i<ntmp;++i) {
        for (size_t j=0;j<refres.size();++j) {
            if (dndt[i]>0) oq<<sqrt(resprops[i][j]/dndt[i])<<"  ";
            else oq<<0.0<<"  ";
        }

        oq<<"\n";
    }

    oq.close();
}

