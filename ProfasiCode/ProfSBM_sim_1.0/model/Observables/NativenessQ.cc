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

#include "NativenessQ.hh"
#include "../Aux/Constants.hh"
#include "../Aux/RMSD_Utils.hh"

using namespace UnivConstants;

using namespace prf;

using namespace prf_utils;
using std::string;
using std::valarray;
using std::vector;
using std::list;

NativenessQ::NativenessQ()
{
    Name("NativenessQ");
    fixed_his=true;
    altcrd="first";
}

NativenessQ::~NativenessQ() {}

/**
\page opt_NativenessQ Nativeness Q parameter
\section options Available Options
\li \b structure structure abcd.pdb:5:A,28,59<br> Reference structure is in
the PDB file abcd.pdb, model 5, chain A, residues index labels 28 through
59 in that file. The NativenessQ %Observable obeys the same rules for
selecting parts of PDB files as prf::ProteinRMSD. The rules are described
in \ref prf_sel_fils .
\li \b alt_crd: alt_crd occupancy \n In case an alternative atom entries
    are found in a pdb file, this command sets the policy to be adopted. The
    value "occupancy" ensures the line with the highest occupancy value is
    used. Other possible policies are "first" and "last".
\sa prf::NativenessQ
*/
int NativenessQ::init_obs()
{
    if (Observable::init_obs()==0) return 0;

    vector<string> parts;

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);

        if ((parts[0]==string("struc") || parts[0]==string("structure"))
            && parts.size()>=2) {
            nativeStructure(parts[1]);
        }
        if ((parts[0]==string("alt_crd"))) altcrd=parts[1];
}

    return init();
}

void NativenessQ::his_range(double xmn, double xmx)
{
    xmin0=0; xmax0=1; //no matter what anyone tries to set these to.
    xbin0=0.01;
}

int NativenessQ::init()
{
    Logger blog(log_thres);
    dcut=6.0;
    d0=3.0;
    scale2=d0*d0;
    blog<<"Initialization of nativeness observable "<<Name()<<"\n";

    if (p==NULL) {
        prf::cerr<<Name()<<": No population set. Initialization failed. \n";
        return 0;
    }

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
        <<"NativenessQ initialization failed.\n";
        return 0;
    }

    list<SelRes> slist1,slist2;

    pdb.mk_selection(sel1,slist1);
    p->mk_selection(sel2,slist2);
    //Sequence alignment, residue listing.
    RMSD_Utils utils;

    if (!utils.seq_match(slist1,slist2)) {
        prf::cerr<<Name()<<": Sequences of the two specified objects "
        <<"do not match.\n Performing sequence alignment...\n";
        utils.seq_align(slist1,slist2);

        if (slist1.empty()||slist2.empty()) {
            prf::cerr<<Name()<<"Attempt to align sequences produced empty"
            <<" lists. "<<Name()<<" initialization failed.\n";
            return 0;
        }
    }

    list<AtomDescriptor> des1,des2;

    pdb.descriptors(slist1,des1);
    p->descriptors(slist2,des2);
    utils.make_filter(string("+CA"));
    utils.apply_filter(des1);
    utils.apply_filter(des2);
    utils.alt_coord_policy(altcrd);
    utils.alt_coord_parse(des1);
    utils.alt_coord_parse(des2);
    nres=des2.size();
    res.resize(nres);
    Qlochist.allocate(ntmp,n_loc_vals());
    Qlochist*=0;
    dndt.resize(ntmp);
    dndt*=0;
    int js2=0;

    for (list<SelRes>::iterator i=slist2.begin();i!=slist2.end();++i) {
        res[js2++]=i->indx_nat;
    }

    //Allocate arrays res, Qloc, n_neighbours, neighbour,refdist
    Qloc.resize(nres,0.0);

    n_neighbours.resize(nres,0);

    neighbour.resize(nres);

    refdist.resize(nres);

    //Loop through res and fill n_neighbours, neighbour
    Shape natv;

    vector<int> vct(des1.size(),0);

    cas.resize(des2.size(),0);

    js2=0;

    for (list<AtomDescriptor>::iterator i=des1.begin();i!=des1.end();++i)  {
        vct[js2++]=i->int_label;
    }

    js2=0;

    for (list<AtomDescriptor>::iterator i=des2.begin();i!=des2.end();++i)  {
        cas[js2++]=i->int_label;
    }

    if (nres!=(int) vct.size() || nres!=(int) cas.size()) {
        prf::cerr<<Name()<<"> Size error ! "
        <<"Native: "<<vct.size()<<", Comparison: "<<cas.size()<<"\n";
        return 0;
    }

    pdb.export_shape(vct,natv);

    for (int i=0;i<nres;++i) {
        n_neighbours[i]=0;
        Vector3 v1=natv.Point(i);

        for (int j=0;j<nres;++j) {
            if (abs(i-j)<3) continue;

            Vector3 v2=natv.Point(j);

            double tmpdist=(v1-v2).mag();

            if (tmpdist<dcut) ++n_neighbours[i];
        }

        if (n_neighbours[i]>0) {
            blog<<"Number of neighbours for residue "<<i<<" is "
            <<n_neighbours[i]<<"\n";
            int jngb=0;
            neighbour[i].resize(n_neighbours[i],0);
            refdist[i].resize(n_neighbours[i],0.0);

            for (int j=0;j<nres;++j) {
                if (abs(i-j)<3) continue;

                Vector3 v2=natv.Point(j);

                double tmpdist=(v1-v2).mag();

                if (tmpdist<dcut) {
                    neighbour[i][jngb]=j;
                    refdist[i][jngb]=tmpdist;
                    ++jngb;
                }
            }
        }
    }

    ntotpairs=0;

    for (int i=0;i<(int) cas.size();++i) {
        for (int j=0;j<n_neighbours[i];++j) {
            if (neighbour[i][j]>i) ++ntotpairs;
        }
    }

    blog<<Name()<<"> Total number of pairs of neighbours = "

    <<ntotpairs<<"\n";

    if (ntotpairs==0) {
        prf::cerr<<Name()<<"> Either there is nothing to do, or "
        <<"the initialization has somehow failed.\n";
        return 0;
    }

    blog<<Name()<<"> Successful initialization with "<<infile<<"\n";

    return 1;
}

double NativenessQ::evaluate()
{
    double obvl=0;
    Shape curshp;
    p->export_shape(cas,curshp);

    for (int i=0;i<curshp.NPoints();++i) {
        Qloc[i]=0;
        Vector3 v1=curshp.Point(i);

        for (int j=0;j<n_neighbours[i];++j) {
            double tmp=(curshp.Point(neighbour[i][j])-v1).mag();
            double qtmp=exp(-(tmp-refdist[i][j])*(tmp-refdist[i][j])/scale2);
            Qloc[i]+=qtmp;

            if (neighbour[i][j]>i) obvl+=qtmp;
        }

        if (n_neighbours[i]!=0) Qloc[i]/=n_neighbours[i];
    }

    return obvl/ntotpairs;
}

void NativenessQ::loc_val(size_t i, int & resi, double &qi)
{
    if (i>0 && i< Qloc.size()) {
        qi=Qloc[i];
        resi=res[i];
    }
}

void NativenessQ::avg_fill(int i)
{
    Observable::avg_fill(i);
    if (gathstat) {
        for (size_t j=0;j<Qloc.size();++j) Qlochist[i][j]+=Qloc[j];
        dndt[i]++;
    }
}

void NativenessQ::avg_reset()
{
    Observable::avg_reset();
    Qlochist*=0;
    dndt*=0;
}

void NativenessQ::avg_write(Output &op)
{
    Observable::avg_write(op);
    Output oq((oprefx+Name()+string(".profile")).c_str(),"w");

    for (int i=0;i<ntmp;++i) {
        for (size_t j=0;j<Qloc.size();++j) {
            if (dndt[i]>0) oq<<Qlochist[i][j]/dndt[i]<<"  "; else oq<<0.0<<"  ";
        }

        oq<<"\n";
    }

    oq.close();
}

