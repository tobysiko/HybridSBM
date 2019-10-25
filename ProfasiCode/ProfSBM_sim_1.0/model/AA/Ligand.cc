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

#include "Ligand.hh"
#include <sstream>

using namespace prf;
using std::string;

Ligand::Ligand(OneLetterCode cod)
{
    mytype=cod;
    rc=lc=NULL;
    mygrp=&Groups::grp[mytype];
    natms=mygrp->num_atoms();
    nhvatm=mygrp->num_hvatoms();
    isaa=mygrp->isAA();
    iseg=mygrp->isEG();
    nrtdof=mygrp->num_sdof();
    unid=atOffset=rtOffset=na=nd=-999;
    chloc=chser=seqser=0;
}

Ligand::Ligand(const Ligand &glg) : rc(glg.rc), lc(glg.lc), unid(glg.unid),
        chloc(glg.chloc), atOffset(glg.atOffset), nrtdof(glg.nrtdof),
        rtOffset(glg.rtOffset), na(glg.na), nd(glg.nd), natms(glg.natms),
        nhvatm(glg.nhvatm), seqser(glg.seqser), chser(glg.chser), atm(glg.atm),
        node(glg.node),
        isaa(glg.isaa), iseg(glg.iseg), mytype(glg.mytype), mygrp(glg.mygrp) {}

Ligand & Ligand::operator=(const Ligand &glg)
{
    if (this!=&glg) {
        mytype=glg.mytype;
        rc=glg.rc;
        lc=glg.lc;
        mygrp=glg.mygrp;
        natms=glg.natms;
        nhvatm=glg.nhvatm;
        node=glg.node;
        isaa=glg.isaa;
        iseg=glg.iseg;
        nrtdof=glg.nrtdof;
        unid=glg.unid;
        atOffset=glg.atOffset;
        rtOffset=glg.rtOffset;
        na=glg.na;
        nd=glg.nd;
        chloc=glg.chloc;
        chser=glg.chser;
        seqser=glg.seqser;
    }

    return *this;
}

Ligand::~Ligand() {}

void Ligand::Initialize() {}

void Ligand::Reconstruct() {}

void Ligand::BuildConnections() {}

void Ligand::atomOffset(int j)
{
    int prevoffset=atm[0].UniqueId();
    for (int i=1;i<numAtoms();++i)
        if (atm[i].UniqueId()<prevoffset) prevoffset=atm[i].UniqueId();
    atOffset=j;

    for (size_t i=0;i<atm.size();++i)
        atm[i].UniqueId(atm[i].UniqueId()-prevoffset +j);
}

Atom * Ligand::labeled_atom(string alabel)
{
    int i=0;
    return (i=mygrp->index(alabel)) <0?NULL:&atm[i];
}

int Ligand::rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1){return 0;}

int Ligand::rotDof_assign(int il, double mgd) {return 0;}

double Ligand::get_rotDof(int il) {return 0;}

int Ligand::ROTDOF(int il, double mgd, int &a0, int &a1) {return 0;}

int Ligand::ROTDOFr(int il, double mgd) {return 0;}

double Ligand::ADOF(int il) {return 0;}

void Ligand::Donor(int i, Dipole &dp) {}

void Ligand::Acceptor(int i, Dipole &dp) {}

void Ligand::ExportConnections(ConnectionsMatrix &aa) {}

void Ligand::LocPairsatRTdof(int i, std::deque<std::pair<int,int> > & lcp) {}

void Ligand::Allocate()
{
    //to resize with natms and run loop otherwise is deliberate!
    atm.resize(natms);

    for (int i=0;i<mygrp->num_atoms();++i) {
        atm[i].Species(mygrp->species(i));
    }
}

Atom & Ligand::at(string alabel)
{
    int indx=mygrp->index(alabel);

    if (indx>=0) return atm[indx];

    prf::cerr<<"Fatal: Atom & Ligand::at(alabel) \n"
    <<"Could not locate atom for label \""<<alabel
    <<"\" on "<<Name() <<"\nAbort.\n";

    exit(1);

    return atm[0];
}

string Ligand::label_of(int i)
{
    return mygrp->label(i);
}

void Ligand::imprint(int iat,AtomDescriptor &ds)
{
    ds.keyword="ATOM";
    ds.iatom=ds.int_label=atom(iat).UniqueId();
    ds.resnm=TLC();
    ds.atom_label=label_of(iat);
    ds.atom_type=ds.atom_label[1];
}

void Ligand::WritePDBline(int &atindx, char ch_id, int aaindx, int i,
                          FILE *fp, int het)
{
    int iat=atm[i].UniqueId();
    AtomCoordinates loc(iat);
    char keywd[][7]={"ATOM  ","HETATM"};
    fprintf(fp,
            "%s%5u %s %s %c%4u    %8.3f%8.3f%8.3f  1.00  0.00           %c  \n",
            keywd[het?1:0],
            ++atindx, (label_of(i)).c_str(), (TLC().c_str()),ch_id,
            aaindx,loc.x(),loc.y(),loc.z(),label_of(i)[1]);

}

void Ligand::WritePDB(int &atindx, char ch_id, int aaindx, FILE *fp){}

void Ligand::WritePDB2(int &atindx, char ch_id, int aaindx, FILE *fp){}

int Ligand::read_coordinates(std::list<AtomRecord>::iterator start,
                             std::list<AtomRecord>::iterator end,
                             std::vector<bool> &assignments)
{
    Atom *a;

    for (; start!=end;++start) {
        if (start->descriptor().keyword!=string("ATOM  ")&&
            start->descriptor().keyword!=string("HETATM")) continue;

        std::string wanted=start->descriptor().atom_label;

        if ((a=labeled_atom(wanted))!=NULL) {
            a->Pos(start->coordinates());
            assignments[a->UniqueId()]=true;
        } else {
            prf::cerr<<"While trying to import record {"<<start->pdb_line()
            <<"} into group "<<Name()<<", could not locate atom \""<<wanted<<"\"!\n";
        }
    }

    return 0;
}

int Ligand::n_dof() const
{
    return n_rotDof();
}

int Ligand::n_coord() const
{
    return n_rotDof();
}

double Ligand::get_coord(int i)
{
    return get_rotDof(i);
}

void Ligand::set_coord(int i, double x)
{
    rotDof_assign(i,x);
}

double Ligand::get_dof(int i)
{
    return get_rotDof(i);
}

void Ligand::set_dof(int i, double x)
{
    rotDof_assign(i,x);
}

int Ligand::set_coord_xml(prf_xml::XML_Node *res, int typecheck)
{
    if (res==NULL or res->name()!="group") return 0;
    if (res->attribute("type")!=TLC() and typecheck) return 0;
    prf_xml::XML_Node *dofs=res->child("coordinates");
    if (!dofs) dofs=res->child("dof");
    if (dofs) {
        std::string values=dofs->value();
        std::istringstream ssin(values);
        double tmpq=0;
        int i=0;
        while (ssin>>tmpq and i<n_coord()) set_coord(i++, tmpq);
    }
    return 1;
}

prf_xml::XML_Node * Ligand::make_xml_node()
{
    prf_xml::XML_Node * xnd=new prf_xml::XML_Node("group","");
    xnd->set_attribute("type",TLC());
    std::ostringstream ost;
    ost<<SeqSerial();
    xnd->set_attribute("index",ost.str());
    ost.str("");
    ost.precision(16);
    for (int i=0;i<n_coord();++i) {
        ost<<get_coord(i);
        if ((i+1)%3 ==0) ost<<"\n"; else ost<<"  ";
    }
    ost<<"\n";
    xnd->add_child_node("coordinates",ost.str());
    return xnd;
}

void Ligand::Write_XML(FILE *fp)
{
    fprintf(fp,"<group index=\"%d\" type=\"%s\">",SeqSerial(), TLC().c_str());
    fprintf(fp,"\n<coordinates>\n");
    for (int i=0;i<n_coord();++i) {
        fprintf(fp,"%.16f %s",get_coord(i),((i+1)%3 ==0)?"\n":"  ");
    }
    fprintf(fp,"\n</coordinates>\n</group>\n");
}

int Ligand::calc_torsions(std::vector<bool> &specified)
{
    int nerr=0;
    int dfatms[4]={0,0,0,0};
    for (int i=0;i<nrtdof;++i) {
        mygrp->dof_def_atoms(i,dfatms);
        int uid[4];
        Vector3 v[4];
        bool specok=true;
        for (int j=0;j<4;++j) {
            if (dfatms[j]<natms) uid[j]=atm[dfatms[j]].UniqueId();
            else uid[j]=dfatms[j]+atOffset;
            v[j]=AtomCoordinates::vec(uid[j]);
            specok=specok&&specified[uid[j]];
        }
        if (specok) {
            Vector3 v1=v[1]-v[0],v2=v[2]-v[1],v3=v[3]-v[2];
            double phai=v2.torsion(v1,v3);
            rotDof_assign(i,phai);
        } else ++nerr;
    }
    return nerr;
}

bool Ligand::get_rotDofAxis(int i, Atom &a0, Atom &a1)
{
    if (i>=0 and i<nrtdof) {
        int dfatms[4]={0,0,0,0};
        mygrp->dof_def_atoms(i,dfatms);
        a0=atm[dfatms[1]];
        a1=atm[dfatms[2]];
        return true;
    }
    return false;
}

Node * Ligand::node_for_dof(int i)
{
    if (i>=0 and i<nrtdof) return node[i];
    else return NULL;
}
