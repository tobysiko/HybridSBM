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

#include "Fragmenter.hh"
#include <sstream>

Fragmenter::Fragmenter() : batchsize(8) {}
Fragmenter::~Fragmenter() {}

void Fragmenter::make_fragments(std::deque<prf_xml::XML_Node *> &segs)
{
    segs.clear();
    for (int i=0;i<p->NumberOfChains();++i) {
        std::deque<prf_xml::XML_Node *> parts, grps;
        Protein *c=p->Chain(i);
        std::ostringstream ost;
        ost<<i;
        std::string origchn=ost.str();
        for (int j=0;j<c->numLigands();++j) {
            prf_xml::XML_Node *grpnd=c->memberLigand(j)->make_xml_node();
            grps.push_back(grpnd);
        }
        size_t iprt=0, ir0=0,ir1=0;
        while (ir1<grps.size()) {
            prf_xml::XML_Node *prtnd=new prf_xml::XML_Node("protein","");
            prtnd->set_attribute("id","0");
            if ((grps.size()-ir0)<=1.5*batchsize) ir1=grps.size();
            else ir1=ir0+batchsize;
            prtnd->set_attribute("original_chain",origchn);
            prtnd->add_child_node("residue_uid_offset",
                                  grps[ir0]->attribute("index"));
            ost.str("");
            ost<<iprt;
            prtnd->add_child_node("segment",ost.str());

            int nres=0;
            if (iprt!=0) {
                prtnd->add_child_node("start_res_orig","1");
                prf_xml::XML_Node *vnd=new prf_xml::XML_Node("group","");
                vnd->set_attribute("type","VoidEG");
                prtnd->add_child_node(vnd);
                ++nres;
            } else prtnd->add_child_node("start_res_orig","0");
            for (size_t ilg=ir0;ilg<ir1;++ilg) {
                prf_xml::XML_Node *vnd=new prf_xml::XML_Node("group","");
                vnd->make_clone_of(grps[ilg]);
                prtnd->add_child_node(vnd);
                ++nres;
            }
            ost.str("");
            ost<<nres;
            prtnd->add_child_node("end_res_orig",ost.str());
            if (!grps.empty()) {
                prf_xml::XML_Node *vnd=new prf_xml::XML_Node("group","");
                vnd->set_attribute("type","VoidEG");
                prtnd->add_child_node(vnd);
            }
            size_t ilg=0;
            std::string seqstr;
            for (size_t ilg1=0;ilg1<prtnd->n_children();++ilg1) {
                if (prtnd->child(ilg1)->name()!="group") continue;
                ost.str("");
                ost<<ilg++;
                prtnd->child(ilg1)->set_attribute("index",ost.str());
                seqstr+=(" "+prtnd->child(ilg1)->attribute("type")+" ");
            }
            prtnd->add_child_node("sequence",seqstr);
            parts.push_back(prtnd);
            ++iprt;
            ir0=ir1-1;
        }
        for (size_t ipr=0;ipr<parts.size();++ipr) {
            prf_xml::XML_Node *tmppop=new prf_xml::XML_Node("population","");
            tmppop->add_child_node("num_chains","1");
            tmppop->add_child_node(parts[ipr]);
            segs.push_back(tmppop);
        }
    }
}

void Fragmenter::join_fragments(std::deque<prf_xml::XML_Node *> &segs)
{
    for (size_t i=0;i<segs.size();++i) {
        prf_xml::XML_Node *prtnd=segs[i]->child("protein");
        int chid=atoi(prtnd->attribute("original_chain").c_str());
        int offset=atoi(prtnd->child("residue_uid_offset")->value().c_str());
        int ires0=atoi(prtnd->child("start_res_orig")->value().c_str());
        int ires1=atoi(prtnd->child("end_res_orig")->value().c_str());
        int ires=0;
        for (size_t j=0;j<prtnd->n_children();++j) {
            prf_xml::XML_Node *lgnd=prtnd->child(j);
            if (lgnd->name()!="group") continue;
            if (ires<ires0 or ires>=ires1) {
                ++ires;
                continue;
            }
            Ligand *lg=p->Chain(chid)->memberLigand(offset+ires-ires0);
            if (lgnd->attribute("type")==lg->TLC()) {
                prf_xml::XML_Node *dofs=lgnd->child("coordinates");
                if (!dofs) dofs=lgnd->child("dof");
                if (dofs) {
                    std::string values=dofs->value();
                    std::istringstream ssin(values);
                    double tmpq=0;
                    std::deque<double> tmpcrds;
                    while (ssin>>tmpq) tmpcrds.push_back(tmpq);
                    int strt=0;
                    if (i!=0 && ires==ires0) strt=1; else strt=0;
                    for (int icrd=strt;icrd<lg->n_coord();++icrd)
                        lg->set_coord(icrd, tmpcrds[icrd]);
                }
            }
            ++ires;
        }
    }
}
