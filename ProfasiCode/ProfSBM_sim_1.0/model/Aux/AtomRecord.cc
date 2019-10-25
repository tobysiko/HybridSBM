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

#include "AtomRecord.hh"
#include "profasi_io.hh"
#include <cstdlib>
#include <map>

using std::string;

namespace prf
{

    using namespace prf_pdb_vars;
    AtomDescriptor::AtomDescriptor()
    {
        int_label=imdl=iatom=0;
        occ="1.00";
        ich=" ";ialt=' ';atom_type='H';
        keyword="KKKKKK";
        resnm="RRR";
        ires="rrr";
        iresrel=0;
        atom_label="LLLL";
    }

    AtomDescriptor::~AtomDescriptor() {}

    AtomDescriptor::AtomDescriptor(const AtomDescriptor &ds)
    {
        int_label=ds.int_label;
        imdl=ds.imdl;
        ires=ds.ires;
        iresrel=ds.iresrel;
        iatom=ds.iatom;
        occ=ds.occ;
        ich=ds.ich;
        ialt=ds.ialt;
        atom_type=ds.atom_type;
        keyword=ds.keyword;
        resnm=ds.resnm;
        atom_label=ds.atom_label;
    }

    AtomDescriptor &AtomDescriptor::operator=(const AtomDescriptor &ds)
    {
        if (this!=&ds) {
            int_label=ds.int_label;
            imdl=ds.imdl;
            ires=ds.ires;
            iresrel=ds.iresrel;
            iatom=ds.iatom;
            occ=ds.occ;
            ich=ds.ich;
            ialt=ds.ialt;
            atom_type=ds.atom_type;
            keyword=ds.keyword;
            resnm=ds.resnm;
            atom_label=ds.atom_label;
        }

        return *this;
    }

    void AtomDescriptor::print_info()
    {
        prf::cout<<int_label<<" "<<" :"<<ires<<":"<<iresrel<<": "<<iatom<<" :"
        <<ich<<": alt:"<<ialt<<": with p "<<occ<<" type "
        <<atom_type<<" "<<resnm<<" "<<atom_label<<"\n";
    }

    string AtomDescriptor::short_info()
    {
        string ans="",sepa="/";
        ans+=ich+sepa+ires+sepa+resnm+sepa+atom_label+sepa;
        return ans;
    }

    void AtomRecord::print_info()
    {
        descriptor().print_info();
        prf::cout<<coordinates()<<"\n";
    }

    bool AtomDescriptor::corresponds_to(AtomDescriptor &ds, string optns)
    {
        return (
                   (atom_type==ds.atom_type) &&
                   (atom_label==ds.atom_label) &&
                   (optns.find("ignore_res_index")<optns.size() || ires==ds.ires) &&
                   (optns.find("ignore_res_label")<optns.size() || resnm==ds.resnm) &&
                   (optns.find("ignore_chain_label")<optns.size() || ich==ds.ich) &&
                   (optns.find("ignore_model_label")<optns.size() || imdl==ds.imdl)
               );
    }

    AtomRecord::AtomRecord() : coord(0.0,0.0,0.0) {}

    AtomRecord::~AtomRecord() {}

    AtomRecord::AtomRecord(string gline, AtomLabelDictionary *dict)
    {
        AtomRecord();
        set_fields(gline,dict);
    }

    /*
    void AtomDescriptor::filter(std::string &alabel)
    {
     std::string dummy="";
     std::map<string,string> lblmap;
     lblmap[" HB2"]="1HB ";
     lblmap[" HB3"]="2HB ";
     lblmap[" HG2"]="1HG ";
     lblmap[" HG3"]="2HG ";
     lblmap[" HA2"]="1HA ";
     lblmap[" HA3"]="2HA ";
        if (alabel[0]=='H') {
            dummy=alabel;

            for (int i=0;i<atom_label_width;++i) {
                dummy[(i+1)%atom_label_width]=alabel[i];
            }

        } else if (lblmap.find(alabel)!=lblmap.end()) {
         dummy=lblmap[alabel];
        }
        if (!dummy.empty()) {
         Logger(20)<<"Interpreting "<<alabel<<" as "<<dummy<<"\n";
         alabel=dummy;
        }
    }
    */

    int AtomDescriptor::set_fields(string gline, AtomLabelDictionary *dict)
    {
        std::string dummy;

        keyword.assign(gline,0,6);

        if (keyword!="ATOM  " && keyword!="HETATM") return 0;

        resnm.assign(gline,res_label_column,res_label_width);

        atom_label.assign(gline,atom_label_column,atom_label_width);

        dict->interpret(resnm,atom_label);
        
        dummy.assign(gline,atom_index_column,atom_index_width);

        iatom=atoi(dummy.c_str());
        atom_type=gline[atom_type_column];
        ialt=gline[alt_coord_column];
        ich.assign(gline,chain_label_column,chain_label_width);
        ires.assign(gline,res_num_column,res_num_width);
        occ=gline.substr(occup_col,occup_width);
        return 1;
    }

    void AtomDescriptor::mark_fields(char gline[81])
    {
        sprintf(gline,"%s%5d %s%c%s %s %s",keyword.c_str(),
                iatom,atom_label.c_str(),ialt,resnm.c_str(),
                ich.c_str(),ires.c_str());

        for (int i=0;i<occup_width;++i) gline[occup_col+i]=occ[i];
    }

    void AtomRecord::build_pdb_line_from_fields()
    {
        char tmpline[81],crd[25];

        for (int i=0;i<81;++i) tmpline[i]=' ';

        dsc.mark_fields(tmpline);

        for (int i=0;i<80;++i) if (tmpline[i]==0) tmpline[i]=' ';

        tmpline[80]=0;

        sprintf(crd,"%8.3f%8.3f%8.3f",coord.x(),coord.y(),coord.z());

        for (int i=0;i<24;++i) {
            tmpline[30+i]=crd[i];
        }

        tmpline[77]=dsc.atom_type;

        char tempfact[]="1.00  0.00";

        for (int i=0;i<10;++i) tmpline[56+i]=tempfact[i];

        line=string(tmpline);
    }

    int AtomRecord::set_fields(string gline, AtomLabelDictionary *dict)
    {
        line=gline;

        if (dsc.set_fields(gline,dict)==0) return 0;

        coord=Vector3(atof(gline.substr(crd_x_col,crd_col_width).c_str()),
                      atof(gline.substr(crd_y_col,crd_col_width).c_str()),
                      atof(gline.substr(crd_z_col,crd_col_width).c_str()));

        return 1;
    }
}
