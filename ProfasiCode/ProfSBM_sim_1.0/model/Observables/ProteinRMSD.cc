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

#include "ProteinRMSD.hh"

using namespace prf_utils;
using std::string;
using std::valarray;
using std::vector;
using std::list;

namespace prf
{
    ProteinRMSD::ProteinRMSD()
    {
        rnge=0;
        fixed1=true;
        calcQ=false;
        seq_align=true;
        altcrd="first";
        dsel2=dsel1="1:*";
        utils.set_logger_threshold(log_thres);
        grdtyp=2;
    }

    ProteinRMSD::~ProteinRMSD() {}

    void ProteinRMSD::rangeEstimate(double &x1, double &x2)
    {
        if (calcQ) {
            x1=0;x2=1;
            if (!userbinsz) xbin0=0.01;
        } else {
            x1=0;x2=rnge;
            if (!userbinsz) xbin0=0.25;
        }
    }

    void ProteinRMSD::filter(string flt)
    {
        myfilter=flt;
        utils.make_filter(myfilter);
    }

    void ProteinRMSD::set_logger_threshold(int thr)
    {
        Observable::set_logger_threshold(thr);
        utils.set_logger_threshold(thr);
    }
    /**
    \page opt_ProteinRMSD ProteinRMSD : RMSD of arbitrary type.
    \section options Available Options
    \li \b struc1: struc1 abcd.pdb:5:A,28,59 \n Reference structure taken from
    the PDB file abcd.pdb, model 5, chain A, residues index labels 28 through 59
    in that file. The file name abcd.pdb can be replaced by the \$ sign, to mean
    that the comparison structure is a live state of the Population in the
.
    \li \b struc2: struc2 efgh.pdb:5:A,28,59 \n The second structure for the
    comparison. The filename can be replaced by the \$ sign, to mean that the
    comparison structure is a live state of the Population in the simulation.
    \li \b using: using +\@C-BB \n Sets the filter to choose atoms for
    RMSD evaluations based on their properties. See \ref prf_sel_fils for
    details.
    \li \b no_align: no_align; A one word command that suppresses any attempt
    to align the sequences. The user is then responsible for all about lining
    up the two sequences.
    \li \b alt_crd: alt_crd occupancy \n In case an alternative atom entries
    are found in a pdb file, this command sets the policy to be adopted. The
    value "occupancy" ensures the line with the highest occupancy value is
    used. Other possible policies are "first" and "last".
    \li \b gradient_method: gradient_method 0 \n To force the use of the slower
    trivial method to calculate gradients: evaluate each partial derivative
    by changing a DOF, reconstructing the population and calculating RMSD
    over and over again. The default method uses Jacobians and is faster. This
    option is only useful for debugging. Any value other than 2 is interpreted
    as a choice to use the slower gradient.
    \sa prf::ProteinRMSD
    */
    int ProteinRMSD::init_obs()
    {
        if (Observable::init_obs()==0) return 0;

        vector<string> parts;

        for (size_t i =0;i<usrcmd.size();++i) {
            parts.clear();
            split(usrcmd[i],parts);

            if (parts[0]==string("struc1") && parts.size()>=2) {
                set_struc1(parts[1]);
            }

            if (parts[0]==string("struc2") && parts.size()>=2) {
                set_struc2(parts[1]);
            }

            if ((parts[0]==string("using")||parts[0]==string("filter"))
                && parts.size()>=2) {
                filter(parts[1]);
            }

            if ((parts[0]==string("no_align"))) no_seq_align();

            if ((parts[0]==string("alt_crd"))) altcrd=parts[1];

        }

        return init();
    }

    int ProteinRMSD::init()
    {
        Logger blog(log_thres);
        string fil1,fil2;
        size_t icolon=fl1.find(':');

        if (icolon<fl1.size()-1) sel1=string(fl1,icolon+1);
        else sel1=dsel1;

        fil1=string(fl1,0,icolon);

        icolon=fl2.find(':');

        if (icolon<fl2.size()-1) sel2=string(fl2,icolon+1);
        else sel2=dsel2;

        fil2=string(fl2,0,icolon);

        if ((fil1=="$" || fil2=="$")&&p==NULL) {
            prf::cerr<<"ProteinRMSD: Evaluation of RMSD with one or both \n"
            <<"structures from the population can only succeed if\n"
            <<"a population object is known to ProteinRMSD.\n"
            <<"ProteinRMSD initialization failed.\n";
            return 0;
        }

        if (fil1=="$") obj1=&(*p);
        else {
            pdb1.set_file(fil1);

            if (pdb1.read_matrix()==0) {
                prf::cerr<<"ProteinRMSD> Could not read in information from "
                <<fl1<<"\n"<<"ProteinRMSD initialization failed.\n";
                return 0;
            }

            obj1=&pdb1;
        }

        if (fil2=="$") obj2=&(*p);
        else {
            pdb2.set_file(fil2);

            if (pdb2.read_matrix()==0) {
                prf::cerr<<"ProteinRMSD> Could not read in information from "
                <<fl2<<"\n"<<"ProteinRMSD initialization failed.\n";
                return 0;
            }

            obj2=&pdb2;
        }

        list<SelRes> slist1,slist2;

        obj1->mk_selection(sel1,slist1);
        obj2->mk_selection(sel2,slist2);

        if (slist1.size()!=slist2.size() && !seq_align) {
            prf::cerr<<"ProteinRMSD> Sequence alignment explicitly disabled. "
            <<"But the two selections have different sizes : "<<slist1.size()
            <<" and "<<slist2.size()<<". Seems like incorrect selections\n"
            <<"Refusing to proceed further with initialisation.\n";
            return 0;
        }

        if (seq_align && !utils.seq_match(slist1,slist2)) {
            blog<<"ProteinRMSD> Sequences of the two specified objects "
            <<"do not match. Performing sequence alignment...\n";
            utils.seq_align(slist1,slist2);

            if (slist1.empty()||slist2.empty()) {
                prf::cerr<<"ProteinRMSD> Attempt to align sequences produced"
                <<" empty lists. ProteinRMSD initialization failed.\n";
                return 0;
            }
        }


        rnge=slist1.size();

        list<AtomDescriptor> des1,des2,mdes1,mdes2;
        obj1->descriptors(slist1,des1);
        obj2->descriptors(slist2,des2);
        utils.apply_filter(des1);
        utils.apply_filter(des2);
        utils.alt_coord_policy(altcrd);
        utils.alt_coord_parse(des1);
        utils.alt_coord_parse(des2);
        make_mapped_lists(des1,des2,mdes1,mdes2);
        sel_atoms1.resize(mdes1.size(),0);
        sel_atoms2.resize(mdes2.size(),0);
        int k=0;

        for (list<AtomDescriptor>::iterator i=mdes1.begin(),j=mdes2.begin();
             i!=mdes1.end();++i,++j,++k) {
            sel_atoms1[k]=i->int_label;
            sel_atoms2[k]=j->int_label;
        }

        if (fixed1) {
            obj1->export_shape(sel_atoms1,shape1);

            if (obj1==&pdb1) pdb1.delete_matrix();

            shape1.find_and_move_to_cm();
        }

        blog<<"ProteinRMSD ("<<Name()<<") with filter "
        <<myfilter<<" between structures "<<fl1<<" and "<<fl2<<"\n";

        if (!seq_align) blog<<Name()
            <<"> No attempt to align sequence according to user's request.\n";

        return 1;
    }

    double ProteinRMSD::evaluate()
    {
        if (!fixed1) {
            obj1->export_shape(sel_atoms1,shape1);
        }

        obj2->export_shape(sel_atoms2,shape2);

        double obvl=rmsd(shape1,shape2,fixed1);

        if (calcQ) return exp(-obvl*obvl/100);
        else return obvl;
    }

    double ProteinRMSD::gradientXYZ(std::valarray<double> &g)
    {
        if (!fixed1) {
            obj1->export_shape(sel_atoms1,shape1);
        }

        obj2->export_shape(sel_atoms2,shape2);

        std::valarray<double> tmpgrd1(3*shape2.NPoints());
        obsval=rmsd.eval_w_grad(shape1,shape2,tmpgrd1);
        g=0;
        for (size_t i=0,j=0;i<sel_atoms2.size();++i) {
            j=sel_atoms2[i];
            g[3*j]=tmpgrd1[3*i];
            g[3*j+1]=tmpgrd1[3*i+1];
            g[3*j+2]=tmpgrd1[3*i+2];
        }
        return obsval;
    }

    void ProteinRMSD::make_mapped_lists(list<AtomDescriptor> &l1,
                                        list<AtomDescriptor> &l2,
                                        list<AtomDescriptor> &ml1,
                                        list<AtomDescriptor> &ml2)
    {
        Logger blog(log_thres);
        list<AtomDescriptor> k1,k2,r1,r2;
        list<AtomDescriptor>::iterator it,jt;
        ml1.clear();ml2.clear();

        while (!(l1.empty()||l2.empty())) {
            k1.clear();k2.clear();r1.clear();r2.clear();
            int opres1=l1.begin()->iresrel;
            int opres2=l2.begin()->iresrel;

            it=l1.begin();jt=l2.begin();

            while (it!=l1.end() && it->iresrel==opres1) ++it;

            while (jt!=l2.end() && jt->iresrel==opres2) ++jt;

            k1.splice(k1.end(),l1,l1.begin(),it);

            k2.splice(k2.end(),l2,l2.begin(),jt);

            for (it=k1.begin();it!=k1.end();) {
                bool match_found=false;
                string matchoptns="ignore_chain_label ignore_model_label ";
                matchoptns+=string("ignore_res_index");

                if (!seq_align) matchoptns+=string("ignore_res_label");

                for (jt=k2.begin();jt!=k2.end();++jt) {
                    if (it->corresponds_to(*jt,matchoptns)) {
                        match_found=true;
                        break;
                    }
                }

                if (match_found) {
                    r2.splice(r2.end(),k2,jt);
                    jt=it;
                    ++jt;
                    r1.splice(r1.end(),k1,it);
                    it=jt;
                } else {
                    blog<<"ProteinRMSD> No match found for "
                    <<it->short_info()<<"\nComparison list was \n";

                    for (jt=k2.begin();jt!=k2.end();++jt) {
                        blog<<jt->short_info()<<"\n";
                    }

                    ++it;
                }
            }

            if (r1.size()==r2.size() && !r1.empty()) {
                ml1.splice(ml1.end(),r1);
                ml2.splice(ml2.end(),r2);
            }
        }

        it=ml1.begin();jt=ml2.begin();

        blog<<"ProteinRMSD> Atom mappings...\n";
        blog<<"List 1 has size "<<ml1.size()<<", and list 2 has size "<<ml2.size()<<"\n";

        while (it!=ml1.end()) {
            blog<<"structure 1 ("<<it->short_info()<<") matched to ";
            blog<<"structure 2 ("<<jt->short_info()<<")\n";
            ++it;++jt;
        }

        prf::clog.flush();
    }
}
