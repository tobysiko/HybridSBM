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

#include "Population.hh"
#include "../Aux/PDBReader.hh"
#include "../Aux/SeqBuilder.hh"
#include "../Aux/RMSD_Utils.hh"
#include "../Aux/profasi_version.hh"
#include "../Aux/prf_time.hh"
#include <fstream>
#include <sstream>

using std::ifstream;
using std::string;
using std::vector;
using std::list;

using namespace prf_utils;

using namespace UnivConstants;

using namespace prf::Groups;

namespace prf
{
    Population::Population() : PopBase(), ntotatms(0),initzd(false)
    {
        charged_NTs=charged_CTs=true;
        RandomNumberGenerator(&triv);
    }

    Population::~Population()
    {
        for (int i=0;i<nc;++i) if (chv[i]) delete chv[i];
    }

    void Population::clear()
    {
        for (size_t i=0;i<chv.size();++i) chv[i]->clear();
        chv.clear();
        lgv.clear();
        aav.clear();
        atyp.clear();
        tmpchv.clear();
        cisres.clear();
        ntotatms=longest=shortest=naam=ndof=0;
        initzd=false;
        charged_NTs=charged_CTs=true;
    }

    void Population::RandomNumberGenerator(RandomNumberBase *rn)
    {
        rnd=rn;
    }

    int Population::num_grp(int i) const
    {
        return chv[i]->numLigands();
    }

    Ligand * Population::existing_group(int ires, int ic)
    {
        return chv[ic]->memberLigand(ires);
    }

    string Population::chain_name(int i) const
    {
        string ans="";

        if (i<0) i=0;

        if (i>=nc) i=nc-1;

        ans.push_back('A'+i);

        return ans;
    }

    int Population::index_of_grp(std::string ires,int ich)
    {
        int ans=atoi(ires.c_str())-1;

        if (ans<0 || ans>=num_grp(ich)) return num_grp(ich);

        return ans;
    }

    string Population::grp_name(int ires, int ic)
    {
        return Chain(ic)->memberLigand(ires)->TLC();
    }

    int Population::descriptors(std::list<SelRes> &slc,
                                std::list<AtomDescriptor> &des)
    {
        int ndes=0;
        des.clear();

        for (list<SelRes>::iterator r=slc.begin();r!=slc.end();++r) {
            Ligand *lg=existing_group(r->indx_nat,r->chn);

            for (int iat=0;iat<lg->numAtoms();++iat) {
                AtomDescriptor d;
                lg->imprint(iat,d);
                d.iresrel=r->indx_nat;
                char tmp[8];
                sprintf(tmp,"%3d",r->indx_nat);
                d.ires=tmp;
                d.ich=chain_name(r->chn);
                des.push_back(d);
                ++ndes;
            }
        }

        return ndes;
    }

    int Population::AddProtein(string ntg,string sq, string ctg,int howmany)
    {
        Logger blog;
        blog(9) <<"Queueing up "<<howmany<<" copies of Protein "
        << (ntg.c_str()) <<"--"<< (sq.c_str()) <<"--"<< (ctg.c_str())
        <<"("<<sq.size() <<" amino acids)"<<"\n";

        while (howmany-->0) {
            tmpchv.push_back("< "+ntg+" * "+sq+" * "+ctg+" >");
        }
        initzd=false;
        return howmany;
    }

    void Population::writePDBHeader(FILE *fp, unsigned long itime, int tindex, double entot)
    {
        unsigned int remn=1;
        fprintf(fp,"REMARK%4u PROGRAM : PROFASI VERSION %s\n",remn++,
                profasi_version().c_str());
        fprintf(fp,"REMARK%4u SNAPSHOT: ENERGY = %f time = %lu\n",remn++,
                entot,itime);
        writeSequenceInfo(fp);
        AtomCoordinates::write_pdb_box_info(fp);
    }

    void Population::writeSequenceInfo(FILE *fp)
    {
        for (int i=0;i<nc;++i) {
            unsigned int remn=1;
            int ilig=0,nlig=chv[i]->numLigands();
            char prol[19];
            do {
                char pepid=' ';
                if (nc<=62) {
                    if (i<26) pepid='A'+i;
                    else if (i<36) pepid='0'+(i-26);
                    else pepid='a'+(i-36);
                }
                sprintf(prol,"SEQRES %3d %c %4d ",remn++,pepid,nlig);
                std::string line=prol;
                for (int j=0;j<13;++j,++ilig) {
                    if ((ilig)<nlig)
                        line+=(" "+chv[i]->memberLigand(ilig)->TLC());
                    else break;
                }
                fprintf(fp,"%s\n",line.c_str());
            } while (ilig<nlig);
        }
    }

    void Population::SaveSnapshot(int in_format, std::string flnm,
                                 unsigned long itime, int tindex, double entot)
    {
        if (in_format==0) return;
        FILE * fp=fopen(flnm.c_str(),"w");
        switch (in_format) {
            case 1: {
                    writePDBHeader(fp,itime,tindex,entot);
                    WritePDB(fp);
                    break;
                }
            case 2: {
                    fprintf(fp,"<?xml version=\"1.0\"?>\n<structure>\n");
                    fprintf(fp,"<box_length>%.16f</box_length>\n",AtomCoordinates::boxL());
                    fprintf(fp,"<energy>%.16f</energy>\n",entot);
                    fprintf(fp,"<snapshot_time> %lu </snapshot_time>\n",itime);
                    fprintf(fp,"<profasi_version>%s</profasi_version>\n",
                            profasi_version().c_str());
                    fprintf(fp,"<creation_time>\nUTC %s</creation_time>\n",
                            prf_time().to_UTC().c_str());
                    Write_XML(fp);
                    fprintf(fp,"</structure>\n");
                    break;
                }
            case 3: {
                    WriteConf_text(fp);
                    break;
                }
            case 4: {
                    WriteConf(fp);
                    break;
                }
            default: break;
        };
        fclose(fp);
    }

    void Population::WritePDB(FILE *fp)
    {
        char pepindx='A';
        int atindx=0, rsindx=0;
        if (nc<=62) {
            int i=0;
            for (;i<nc&&i<26;++i) {
                pepindx='A'+i;
                rsindx=0;
                chv[i]->WritePDB(atindx,fp,pepindx,rsindx);
            }
            for (;i<nc&&i<36;++i) {
                pepindx='0'+(i-26);
                rsindx=0;
                chv[i]->WritePDB(atindx,fp,pepindx,rsindx);
            }
            for (;i<nc&&i<62;++i) {
                pepindx='a'+(i-36);
                rsindx=0;
                chv[i]->WritePDB(atindx,fp,pepindx,rsindx);
            }
        } else {
            pepindx=' ';
            for (int i=0;i<nc;++i) {
                chv[i]->WritePDB(atindx,fp,pepindx,rsindx);
            }
        }
    }

    void Population::WritePDB2(FILE *fp)
    {
        char pepindx='A';
        int atindx=0, rsindx=0;
        if (nc<=62) {
            int i=0;
            for (;i<nc&&i<26;++i) {
                pepindx='A'+i;
                rsindx=0;
                chv[i]->WritePDB2(atindx,fp,pepindx,rsindx);
            }
            for (;i<nc&&i<36;++i) {
                pepindx='0'+(i-26);
                rsindx=0;
                chv[i]->WritePDB2(atindx,fp,pepindx,rsindx);
            }
            for (;i<nc&&i<62;++i) {
                pepindx='a'+(i-36);
                rsindx=0;
                chv[i]->WritePDB2(atindx,fp,pepindx,rsindx);
            }
        } else {
            pepindx=' ';
            for (int i=0;i<nc;++i) {
                chv[i]->WritePDB2(atindx,fp,pepindx,rsindx);
            }
        }
    }

    int Population::re_index()
    {
        Logger blog;
        std::vector<Protein *> tm=chv;
        chv.clear();
        for (size_t i=0;i<tm.size();++i) {
            if (tm[i]!=NULL) chv.push_back(tm[i]);
        }
        nc=(int)chv.size();
        ntotatms=0;
        int iaast=0,sz=0,szs=100000,szl=-1,naa_tot=0;
        longest=shortest=0;

        for (int i=0;i<nc;++i) {
            chv[i]->setAtomOffset(ntotatms);
            ntotatms+=chv[i]->numberOfAtoms();
            iaast+= (sz=chv[i]->numLigands());
            naa_tot+=chv[i]->numAminoAcids();
            if (sz>szl) {szl=sz;longest=i;}
            if (sz<szs) {szs=sz;shortest=i;}
        }

        AtomCoordinates::Initialize(ntotatms);
        ndof=0;
        for (int i=0;i<nc;++i) ndof+=chv[i]->n_dof();
        lgv.resize(iaast,NULL);
        aav.resize(naa_tot,NULL);
        int iaac=0,iaac2=0;

        for (int i=0;i<nc;++i) {
            chv[i]->setUniqueId(i);
            for (int j=0;j<chv[i]->numLigands();++j,++iaac) {
                lgv[iaac]=chv[i]->memberLigand(j);
                lgv[iaac]->UniqueId(iaac);
            }

            for (int j=0;j<chv[i]->numAminoAcids();++j,++iaac2) {
                aav[iaac2]=chv[i]->AA(j);
            }
        }

        atyp.resize(ntotatms,hydrogen);

        vector<bool> covered(ntotatms,false);

        for (unsigned int i=0;i<lgv.size();++i) {
            for (int k=0;k<lgv[i]->NumberOfAtoms();++k) {
                Atom a=lgv[i]->atom(k);
                blog(100)<<"Ligand id "<<i<<" atom "<<k<<" has UID "
                <<a.UniqueId()<<" and species "<<a.Species()<<"\n";
                atyp[a.UniqueId()]=a.Species();
                covered[a.UniqueId()]=true;
            }
        }

        bool allcovered=true;

        for (size_t i=0;i<covered.size();++i) {
            if (!covered[i]) {
                allcovered=false;
                prf::cerr<<"Species of atom "<<i<<" not assigned.\n";
            }
        }

        if (!allcovered) {
            prf::cerr<<"Error in Population initialization.\n"
            <<"Atom type array in Population has not been \n"
            <<"correctly assigned for all atoms.\n";
        }

        calcNSpecies();
        blog(10) <<"Population> Contents ...\n";
        for (int i=0;i<nc;++i) {
            blog <<"Chain "<<i<<": "<<chv[i]->Sequence()
                    <<" ("<<chv[i]->numAminoAcids() <<" residues, "
                    <<chv[i]->numberOfAtoms()<<" atoms.)\n";
        }
        blog <<"Total number of atoms = "<<ntotatms<<"\n";

        index_dof();
        dofindex.reset();
        dofindex.init(dof_info);
        if (check_DOF_index()==0) {
            prf::cerr<<"Inconsistent index for degrees of freedom in population.\n"
                    <<"Initialisation failed.\n";
            return 0;
        }

        if (nc>0) {
            initzd=true;
            return 1;
        } else {
            initzd=false;
            return 0;
        }
    }

    int Population::Init()
    {
        if (tmpchv.empty()) return 0;

        for (size_t i=0;i<tmpchv.size();++i) {
            chv.push_back(new Protein());
            chv.back()->charged_end(charged_NTs, charged_CTs);
            if (!cisres[i].empty()) {
                for (size_t j=0; j<cisres[i].size(); ++j) {
                    chv.back()->setCis(cisres[i][j]);
                }
            }

            chv.back()->Allocate(0,0,tmpchv[i]);
            if (!chv[i]->internally_consistent()) {
                delete chv.back();
                chv.back()=NULL;
                prf::cerr<<"Population> Could not make a chain out of "
                        <<tmpchv[i]<<"\n";
            }
        }
        tmpchv.clear();
        return re_index();
    }

    void Population::Initialize(int inittyp)
    {
        if (inittyp==0) {
            Randomize();
        } else {
            for (int i=0;i<nc;++i) {
                chv[i]->stretched_init();
            }

            RandomizeRelConf();
        }
    }

    int Population::InitCoord(string init_type)
    {
        bool coordsok=false;
        Logger blog(10);
        init_type=trim_str(init_type);
        blog<<"Population> InitCoord : ("<<init_type<<")\n";
        if (init_type=="none") coordsok=true;
        else if (init_type.substr(0,7)==string("file://")) {
            std::string configfile=init_type.substr(7);
            if (TestFile_r(configfile)) {
                int lastdot=configfile.find_last_of('.');
                std::string extension=configfile.substr(lastdot+1);
                if (extension=="tconf" or extension=="tcnf") {
                    FILE *fp=fopen(configfile.c_str(),"r");
                    ReadConf_text(fp);
                    fclose(fp);
                    coordsok=true;
                } else if (extension=="xml") {
                    prf_xml::XML_Node *pnode=prf_xml::get_xml_tree(configfile);
                    assign_structures(pnode);
                    coordsok=true;
                    if (pnode) delete pnode;
                }
                blog<<"Population> Initial conformation read from "
                <<configfile<<"\n";
            } else {
                prf::cerr<<"Could not open file "<<configfile
                <<" for reading\n";
            }
        } else if (init_type.find("stretched")<init_type.size()) {
            for (int i=0;i<nc;++i) {
                chv[i]->stretched_init();
            }

            coordsok=true;

            blog<< "Population> Stretched initial conformation for chains.\n";
        }

        if (!coordsok) {
            blog << "Population> Random initial conformation.\n";
            Randomize();
        } else if (init_type.find("random_rel")<init_type.size()) {
            blog << "Population> Randomizing initial relative positions "
            <<"and orientations.\n";
            RandomizeRelConf();
        }
        if (coordsok) return 1; else return 0;
    }

    int Population::index_dof()
    {
        dof_info.clear();
        ndof=0;
        for (int i=0;i<nc;++i) ndof+=chv[i]->n_dof();
        dof_info.resize(ndof);
        int idof=0,idofrgd=0,idofrot=0,idofbb=0;
        for (int i=0;i<nc;++i) {
            int ilocdof=0,ilocdofrgd=0,ilocdofrot=0,ilocdofbb=0;
            for (int j=0;j<9;++j) {
                DOF_Info tmp;
                tmp.dof_kind=rigid_body_xyz;
                tmp.group=Chain(i)->first_AA()->UniqueId();
                tmp.specific_index_in_group=j;
                tmp.chain=i;
                tmp.global_index=idof;
                tmp.index_in_chain=ilocdof++;
                tmp.specific_global_index=idofrgd++;
                tmp.specific_index_in_chain=ilocdofrgd++;
                dof_info[idof++]=tmp;
            }
            int prevgrp=-1,grpidx=0;
            for (int j=0;j<chv[i]->numBBdof();++j) {
                DOF_Info tmp;
                tmp.dof_kind=backbone_torsion_angle;
                tmp.chain=i;
                tmp.group=chv[i]->residue_with_bb_dof(j)->UniqueId();
                if (tmp.group!=prevgrp) {
                    grpidx=0;
                    prevgrp=tmp.group;
                } else grpidx++;
                tmp.specific_index_in_group=grpidx;
                tmp.global_index=idof;
                tmp.index_in_chain=ilocdof++;
                tmp.specific_global_index=idofbb++;
                tmp.specific_index_in_chain=ilocdofbb++;
                dof_info[idof++]=tmp;
            }
            prevgrp=-1;
            grpidx=0;
            for (int j=0;j<chv[i]->numRTdof();++j) {
                DOF_Info tmp;
                tmp.dof_kind=sidechain_torsion_angle;
                tmp.chain=i;
                tmp.group=chv[i]->residue_with_rt_dof(j)->UniqueId();
                if (tmp.group!=prevgrp) {
                    grpidx=0;
                    prevgrp=tmp.group;
                } else grpidx++;
                tmp.specific_index_in_group=grpidx;
                tmp.global_index=idof;
                tmp.index_in_chain=ilocdof++;
                tmp.specific_global_index=idofrot++;
                tmp.specific_index_in_chain=ilocdofrot++;
                dof_info[idof++]=tmp;
            }
        }
        return 0;
    }


    double Population::get_dof(size_t i)
    {
        return get_dof(dof_info[i%ndof]);
    }

    void Population::set_dof(size_t i, double vl)
    {
        set_dof(dof_info[i],vl);
    }

    DOF_Info & Population::get_dof_info(size_t ich, size_t i)
    {
        int id=dofindex.get_uid_from_chain_index(ich,i);
        return dof_info[id];
    }

    double Population::get_dof(DOF_Info &d)
    {
        return Chain(d.chain)->get_dof(d.index_in_chain);
    }

    void Population::get_dof(std::vector<double> &dar)
    {
        if (dar.size() < ((size_t) n_dof())) dar.resize(n_dof(),0);
        for (int i=0;i<n_dof();++i) {
            dar[i]=get_dof(i);
        }
    }

    void Population::set_dof(std::vector<double> &dar)
    {
        size_t minsz=n_dof();
        if (dar.size()<minsz) minsz=dar.size();
        for (size_t i=0;i<minsz;++i) set_dof(i,dar[i]);
    }

    void Population::set_dof(DOF_Info &d, double vl)
    {
        Chain(d.chain)->set_dof(d.index_in_chain,vl);
    }

    void Population::set_dof(std::string dofstr, double vl)
    {
        int uid=get_dof_id(dofstr);
        if (uid>=0 and uid<(int)dof_info.size()) {
            set_dof(uid,vl);
        }
    }

    int Population::get_dof_id(std::string dofstr)
    {
        //Analyze syntax
        int nsp=0,nerr=0,uid=-1;
        for (size_t i=0;i<dofstr.size();++i) {
            if (isspace(dofstr[i])) ++nsp;
        }
        if (nsp!=0) {
            prf::cerr<<"DOF identifiers can not contain spaces.\n";
            ++nerr;
        }
        std::vector<std::string> parts;
        prf_utils::split_str(dofstr,':',parts);
        if (parts.size()!=4) {
            prf::cerr<<"Dof identifiers must have 3 ':' characters.\n";
            ++nerr;
        }
        int dtype=-1,ich=-1,ires=-1,idx=-1;
        if (nerr==0) {
            if (!parts[0].empty()) ich=atoi(parts[0].c_str());
            if (!parts[1].empty()) ires=atoi(parts[1].c_str());
            if (!parts[3].empty()) idx=atoi(parts[3].c_str());
            if (parts[2]=="0" or parts[2]=="r" or parts[2]=="rigid_body_xyz")
                dtype=0;
            else if (parts[2]=="1" or parts[2]=="b" or
                     parts[2]=="backbone_torsion_angle") dtype=1;
            else if (parts[2]=="2" or parts[2]=="s" or
                     parts[2]=="sidechain_torsion_angle") dtype=2;
            if (dtype<0 and (ires>=0)) {
                prf::cerr<<"Can not specify residue without specifying DOF type.\n";
                ++nerr;
            }
        }
        if (nerr!=0) {
            prf::cerr<<"Can not interpret \""<<dofstr
                    <<"\" as a DOF id because of "<<nerr<<" errors.\n";
            return -1;
        }

        if (ich<0) {
            // index is global, residue index is global
            if (ires<0) {
                // global index
                if (dtype<0) uid=idx;
                else uid=dofindex.get_uid_from_type_index(dtype,idx);
            } else {
                // index is in residue and according to type
                uid=dofindex.get_uid_from_residue_type(ires,dtype,idx);
            }
        } else {
            // index is local to chain, residue is local to chain
            if (ires<0) {
                // index relative to chain
                if (dtype<0) uid=dofindex.get_uid_from_chain_index(ich,idx);
                else uid=dofindex.get_uid_from_chain_type(ich,dtype,idx);
            } else {
                // must reinterpret residue inside the chain
                if (ich<NumberOfChains() and
                    ires<((int) Chain(ich)->numLigands())) {
                    uid=dofindex.get_uid_from_residue_type(Chain(ich)->memberLigand(ires)->UniqueId(),dtype,idx);
                }
            }
        }

        return uid;
    }

    void Population::Randomize()
    {
        RandomizeIntConf();
        RandomizeRelConf();
    }

    void Population::RandomizeIntConf()
    {
        RandomizeIntConf(0,nc);
    }

    void Population::RandomizeIntConf(int ich)
    {
        chv[ich]->randomize(rnd);
    }

    void Population::RandomizeIntConf(int ich, int jch)
    {
        for (int i=ich;i<jch;++i) RandomizeIntConf(i);
    }

    void Population::RandomizeRelConf()
    {
        RandomizeRelConf(0,nc);
    }

    void Population::RandomizeRelConf(int ich)
    {
        if (nc<=1) return;

        double trmag=AtomCoordinates::boxL();

        Vector3 rotdirc(rnd->shoot(),rnd->shoot(),rnd->shoot());

        Vector3 rotorg(chv[ich]->CenterOfMass());

        double rotmag=twoPi*rnd->shoot();

        AtomCoordinates::BlockRotate(rotmag,chv[ich]->begin_atom(),
                                     chv[ich]->end_atom(),rotorg,rotdirc);

        Vector3 trv(rnd->shoot()-0.5,rnd->shoot()-0.5,rnd->shoot()-0.5);

        trv*= (trmag);

        //cout <<"Translating peptide by "<<trv<<" after initialization\n";
        AtomCoordinates::BlockTranslate(trv,chv[ich]->begin_atom(),
                                        chv[ich]->end_atom());

        chv[ich]->EnforceBC();

        AtomCoordinates::update(chv[ich]->begin_atom(),chv[ich]->end_atom());
    }

    void Population::RandomizeRelConf(int ich, int jch)
    {
        for (int i=ich;i<jch;++i) RandomizeRelConf(i);
    }

    void Population::Reconstruct()
    {
        for (int i=0;i<nc;++i) chv[i]->reconstruct();
    }

    void Population::Write()
    {
        for (int i=0;i<nc;++i) {
            chv[i]->Write();
        }
    }

    void Population::WriteShort()
    {
        if (NumberOfChains() ==0) {
            prf::cout<<"Empty Population.\n";
            return;
        }

        highlight("Population");

        for (int i=0;i<nc;++i) {
            prf::cout<<"Chain "<<i<<" : "<<chv[i]->Sequence() <<"\n";
        }

        prf::cout<<"Number of Chains = "<<NumberOfChains() <<"\n";
        prf::cout<<"Number of Residues = "<<NumberOfResidues() <<"\n";
        prf::cout<<"Number of Atoms = "<<NumberOfAtoms() <<"\n";
    }

    int Population::calcNSpecies()
    {
        naam=1;

        for (int i=1;i<nc;++i) {
            int match=0;

            for (int j=0;j<i;++j) {
                if (chv[i]->isSameTypeAs(*chv[j]))
                    {match=1;break;}
            }

            if (match==0) ++naam;
        }

        Logger()(20) <<"Population consists of "<<naam
                <<" distinct species of peptides\n";
        return naam;
    }

    int Population::AddProtein(string fullseq, int hwmny)
    {
        Logger blog;

        blog(9) <<"Queuing up "<<hwmny<<" copies of protein "<<fullseq<<"\n";

        while (hwmny--) tmpchv.push_back(fullseq);
        initzd=false;
        return 1;
    }

    int Population::AddProtein(int hwmny, string pdbfilename)
    {
        Logger blog;
        list<SelRes> slct;
        size_t icolon=pdbfilename.find(':');
        string filenm(pdbfilename,0,icolon);
        string selections="";

        if (icolon<(pdbfilename.size()-1))
            selections=pdbfilename.substr(icolon+1);
        else selections="1:A";
        if (!TestFile_r(filenm.c_str())) {
            prf::cerr<<"Requested addition of "<<hwmny<<" copies of "
                <<pdbfilename<<" failed!\n";
            return 0;
        }
        PDBReader pdb(filenm.c_str());

        if (pdb.read_matrix() ==0) {
            prf::cerr<<"Could not read in information from "<<filenm<<"\n";
            return 0;
        }

        pdb.mk_selection(selections,slct);
        initzd=false;
        return AddProtein(slct,hwmny);
    }

    int Population::AddProtein(list<SelRes> &slct, int hwmny)
    {
        Logger blog;
        RMSD_Utils utl;
        int ichprev=-1234567,nadd=0;
        string tmpseq="< ";
        list<SelRes> tmpres,slctnew;
        blog(20) <<"AddProtein from selection list of size "<<slct.size() <<"\n";

        for (list<SelRes>::iterator it=slct.begin();it!=slct.end();++it) {
            if (it->chn!=ichprev) {
                if (ichprev!=-1234567 && (!tmpres.empty())) {
                    if (utl.clean_seq(tmpres) && (!tmpres.empty())) {
                        tmpseq="< ";
                        for (list<SelRes>::iterator lr=tmpres.begin();lr!=tmpres.end();++lr) {
                            tmpseq+=(lr->nam)+" ";
                        }
                        tmpseq+=" >";
                        blog(9) <<"Queueing up "<<hwmny<<" copies of chain "<<tmpseq<<"\n";

                        for (int i=0;i<hwmny;++i) {++nadd;tmpchv.push_back(tmpseq);}
                        slctnew.splice(slctnew.end(),tmpres);
                    }
                }
                tmpres.clear();
                ichprev=it->chn;
            }
            tmpres.push_back(*it);
        }
        if (utl.clean_seq(tmpres)&&(!tmpres.empty())) {
            tmpseq="< ";
            for (list<SelRes>::iterator lr=tmpres.begin();lr!=tmpres.end();++lr) {
                tmpseq+=(lr->nam)+" ";
            }
            tmpseq+=" >";
            blog(9) <<"Queueing up "<<hwmny<<" copies of chain "<<tmpseq<<"\n";

            for (int i=0;i<hwmny;++i) {++nadd;tmpchv.push_back(tmpseq);}
            slctnew.splice(slctnew.end(),tmpres);
        }
        slct=slctnew;
        initzd=false;
        return nadd;
    }

    int Population::ImportStructure(std::list<AtomRecord> &pdbrcd,
                                    std::vector<bool> &assignments, int at_chain)
    {
        list<AtomRecord>::iterator chbeg,chend;

        if (pdbrcd.empty()) return 0;
        Logger blog(3);
        chbeg=chend=pdbrcd.begin();
        blog(20)<<"Population> Importing coordinate information from an atom record "
                <<"list of size "<<pdbrcd.size()<<"\n";

        std::string current_chain=chbeg->descriptor().ich;

        while (chend != pdbrcd.end()) {
            if (at_chain>=nc) {
                blog<<"Population> Aborting coordinate import loop as the "
                       <<"target chain number "<<at_chain<<" is out of range.\n";
                break;
            }
            ++chend;
            if (chend==pdbrcd.end() or chend->descriptor().keyword=="TER   "
                or chend->descriptor().ich!=current_chain) {
                blog(20)<<"Importing coordinates from chain "
                        <<current_chain<<" of the atom record list into"
                        <<" chain "<<at_chain<<" of the population. \n";
                int atlg=Chain(at_chain)->first_ligand()->OLC()==VOIDEG?1:0;
                Chain(at_chain)->importXYZ(chbeg,chend,assignments,atlg);
                blog<<"Finished import "<<current_chain<<" --> "
                       <<at_chain++<<"\n";
                while (chend!=pdbrcd.end() && chend->descriptor().keyword!="ATOM  "
                       &&chend->descriptor().keyword!="HETATM") ++chend;
                chbeg=chend;

                if (chbeg!=pdbrcd.end()) current_chain=chbeg->descriptor().ich;
            }
        }

        return 0;
    }

    int Population::guess_missing_coordinates(std::vector<bool> &assignments)
    {
        for (int i=0;i<NumberOfChains();++i) {
            Chain(i)->guess_missing_coordinates(assignments);
        }

        return 0;
    }

    int Population::export_descriptors(list<AtomDescriptor> &lst)
    {
        Ligand *lg;

        for (int ich=0,res=0;ich<NumberOfChains();++ich) {
            for (int ilg=0;ilg<Chain(ich)->numLigands();++ilg) {
                if (!(lg=ligand(ilg))) continue;
                else ++res;

                for (int iat=0;iat<lg->numAtoms();++iat) {
                    AtomDescriptor d;
                    lg->imprint(iat,d);
                    d.iresrel=res;
                    char tmp[8];
                    sprintf(tmp,"%3d",res);
                    d.ires=tmp;
                    d.ich=chain_name(ich);
                    lst.push_back(d);
                }
            }
        }

        return (int) lst.size();
    }

    int Population::export_shape(vector<int> &vct, Shape &shp)
    {
        shp.Resize(vct.size());

        for (size_t i=0;i<vct.size();++i)
            shp.Point(i,AtomCoordinates::vec(vct[i]));

        return (int) vct.size();
    }

    prf_xml::XML_Node * Population::make_xml_node()
    {
        prf_xml::XML_Node *xnd=new prf_xml::XML_Node("population","");
        std::ostringstream ost;
        ost<<NumberOfChains();
        xnd->add_child_node("num_chains",ost.str());
        for (int i=0;i<NumberOfChains();++i) {
            prf_xml::XML_Node *prtnd=Chain(i)->make_xml_node();
            xnd->add_child_node(prtnd);
        }
        return xnd;
    }

    void Population::Write_XML(FILE *op)
    {
        fprintf(op,"<population>\n");
        fprintf(op,"<num_chains>%d</num_chains>\n",NumberOfChains());
        for (int i=0;i<NumberOfChains();++i) Chain(i)->Write_XML(op);
        fprintf(op,"</population>\n");
    }

    int Population::assign_sequences(prf_xml::XML_Node * pnode)
    {
        if (pnode==NULL or pnode->name()!="population") return 0;
        unsigned int xmlnch =
                (unsigned) std::atoi(pnode->child("num_chains")->value().c_str());
        if (xmlnch>chv.size()) {
            for (size_t i=chv.size();i<xmlnch;++i) chv.push_back(NULL);
        }
        for (size_t i=0;i<pnode->n_children();++i) {
            prf_xml::XML_Node *prt=pnode->child(i);
            if (prt->name()!="protein") continue;
            int id=atoi(prt->attribute("id").c_str());
            if (id<0 or id>= (int) xmlnch) {
                prf::cerr<<"protein id "<<id<<" in XML node inconsistent "
                        <<"with population characteristic num_chains = "
                        <<xmlnch<<"\n";
            } else {
                if (chv[id]==NULL) chv[id]=new Protein();
                chv[id]->metamorphose(prt);
            }
        }
        re_index();
        return 1;
    }

    int Population::assign_structures(prf_xml::XML_Node * pnode)
    {
        if (pnode==NULL or pnode->name()!="population") return 0;
        for (size_t i=0;i<pnode->n_children();++i) {
            prf_xml::XML_Node *prt=pnode->child(i);
            if (prt->name()!="protein") continue;
            int id=atoi(prt->attribute("id").c_str());
            if (id<0 or id>=nc) {
                prf::cerr<<"Rejecting protein id "<<id<<" for dof assignment\n";
            } else {
                chv[id]->assign_dof(prt);
            }
        }
        prf_xml::XML_Node * dofn=pnode->child("dof_assignments");
        if (dofn!=NULL) {
            dofn->interpret_formatted_data();
            for (size_t i=0;i<dofn->n_children();++i) {
                prf_xml::XML_Node *adof=dofn->child(i);
                if (adof->name()=="dof") {
                    std::string dofname=adof->attribute("id");
                    double vl=strtod(adof->value().c_str(),NULL);
                    set_dof(dofname,vl);
                }
            }
        }
        return 1;
    }

    int Population::Read_XML(prf_xml::XML_Node *pnode)
    {
        if (assign_sequences(pnode)) return assign_structures(pnode);
        return 0;
    }

    int Population::check_DOF_index()
    {
        DOF_Info tmp;
        int nerr=0;
        for (size_t i=0;i<chv.size();++i) {
            // chain level index of the dof
            for (int j=0;j<(int)Chain(i)->n_dof();++j) {
                int uid=dofindex.get_uid_from_chain_index(i,j);
                if (uid<0) {
                    prf::cerr<<"Could not retrieve dof "<<j
                            <<" for chain "<<i<<" from index\n";
                    ++nerr;
                } else {
                    tmp=get_dof_info(uid);
                    if (tmp.index_in_chain!=j) {
                        prf::cerr<<"Serial number mismatch with index for "
                                <<" dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.chain!=(int)i) {
                        prf::cerr<<"chain id mismatch with index for rigid "
                                <<"body dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                }
            }
            // chain level index within dofs
            // rigid body
            int nrgdof=9;
            for (int j=0;j<nrgdof;++j) {
                int uid=dofindex.get_uid_from_chain_type(i,rigid_body_xyz,j);
                if (uid<0) {
                    prf::cerr<<"Could not retrieve rigid body dof "<<j
                            <<" for chain "<<i<<" from index\n";
                    ++nerr;
                } else {
                    tmp=get_dof_info(uid);
                    if (tmp.dof_kind!=rigid_body_xyz) {
                        prf::cerr<<"DOF kind mismatch with index for rigid "
                                <<"body dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.chain!=(int)i) {
                        prf::cerr<<"chain id mismatch with index for rigid "
                                <<"body dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.specific_index_in_chain!=j) {
                        prf::cerr<<"Serial number mismatch with index for rigid "
                                <<"body dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                }
            }
            // backbone
            int nbbdof=Chain(i)->numBBdof();
            for (int j=0;j<nbbdof;++j) {
                int uid=dofindex.get_uid_from_chain_type(i,backbone_torsion_angle,j);
                if (uid<0) {
                    prf::cerr<<"Could not retrieve backbone dof "<<j
                            <<" for chain "<<i<<" from index\n";
                    ++nerr;
                } else {
                    tmp=get_dof_info(uid);
                    if (tmp.dof_kind!=backbone_torsion_angle) {
                        prf::cerr<<"DOF kind mismatch with index for "
                                <<"backbone dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.chain!=(int)i) {
                        prf::cerr<<"chain id mismatch with index for "
                                <<"backbone dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.specific_index_in_chain!=j) {
                        prf::cerr<<"Serial number mismatch with index for "
                                <<"backbone dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                }
            }
            //sidechain
            int nscdof=Chain(i)->numRTdof();
            for (int j=0;j<nscdof;++j) {
                int uid=dofindex.get_uid_from_chain_type(i,sidechain_torsion_angle,j);
                if (uid<0) {
                    prf::cerr<<"Could not retrieve sidechain dof "<<j
                            <<" for chain "<<i<<" from index\n";
                    ++nerr;
                } else {
                    tmp=get_dof_info(uid);
                    if (tmp.dof_kind!=sidechain_torsion_angle) {
                        prf::cerr<<"DOF kind mismatch with index for "
                                <<"sidechain dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.chain!=(int)i) {
                        prf::cerr<<"chain id mismatch with index for "
                                <<"sidechain dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                    if (tmp.specific_index_in_chain!=j) {
                        prf::cerr<<"Serial number mismatch with index for "
                                <<"sidechain dof "<<j<<" on chain "<<i<<"\n";
                        ++nerr;
                    }
                }
            }

            //residue level index for side chain angles
            for (int j=0;j<Chain(i)->numLigands();++j) {
                for (int k=0;k<Chain(i)->memberLigand(j)->n_rotDof();++k) {
                    int uid=dofindex.get_uid_from_residue_type(
                            Chain(i)->memberLigand(j)->UniqueId(),
                            sidechain_torsion_angle,k);
                    if (uid<0) {
                        prf::cerr<<"Could not retrieve side chain dof "<<k
                                <<" for ligand id "<<Chain(i)->memberLigand(j)->UniqueId()
                                <<" from index\n";
                        ++nerr;
                    } else {
                        tmp=get_dof_info(uid);
                        if (tmp.dof_kind!=sidechain_torsion_angle) {
                            prf::cerr<<"DOF kind mismatch with index for sidechain "
                                    <<"dof "<<j<<" on ligand id "
                                    <<Chain(i)->memberLigand(j)->UniqueId()<<"\n";
                            ++nerr;
                        }
                        if (tmp.specific_index_in_group!=k) {
                            prf::cerr<<"Serial number mismatch with index for "
                                    <<"sidechain dof "<<k<<" on ligand id "
                                    <<Chain(i)->memberLigand(j)->UniqueId()<<"\n";
                            ++nerr;
                        }
                        if (tmp.group!=Chain(i)->memberLigand(j)->UniqueId()) {
                            prf::cerr<<"ligand id mismatch for group "<<j<<" for "
                                    <<"sidechain dof on chain "<<i<<"\n";
                            ++nerr;
                        }
                    }
                }
            }

            //residue level index for backbone angles
            for (int j=0;j<Chain(i)->numAminoAcids();++j) {
                int nbbdof= (Chain(i)->AA(j)->OLC()==P or
                             Chain(i)->AA(j)->OLC()==DPR)?1:2;
                for (int k=0;k<nbbdof;++k) {
                    int uid=dofindex.get_uid_from_residue_type(
                            Chain(i)->AA(j)->UniqueId(),
                            backbone_torsion_angle,k);
                    if (uid<0) {
                        prf::cerr<<"Could not retrieve backbone dof "<<k
                                <<" for amino acid with uid "<<Chain(i)->AA(j)->UniqueId()
                                <<" from index\n";
                        ++nerr;
                    } else {
                        tmp=get_dof_info(uid);
                        if (tmp.dof_kind!=backbone_torsion_angle) {
                            prf::cerr<<"DOF kind mismatch with index for backbone "
                                    <<"dof "<<j<<" on amino acid uid "
                                    <<Chain(i)->AA(j)->UniqueId()<<"\n";
                            ++nerr;
                        }
                        if (tmp.specific_index_in_group!=k) {
                            prf::cerr<<"Serial number mismatch with index for "
                                    <<" backbone dof "<<k<<" on aminoacid uid "
                                    <<Chain(i)->AA(j)->UniqueId()<<"\n";
                            ++nerr;
                        }
                        if (tmp.group!=Chain(i)->AA(j)->UniqueId()) {
                            prf::cerr<<"ligand id mismatch for group "<<j<<" for "
                                    <<"backbone dof on chain "<<i<<"\n";
                            ++nerr;
                        }
                    }
                }
            }
        }

        // system wide index within a category of DOFs
        //rigid body DOFs
        int nrgdof=9*chv.size();
        for (int j=0;j<nrgdof;++j) {
            int uid=dofindex.get_uid_from_type_index(rigid_body_xyz,j);
            if (uid<0) {
                prf::cerr<<"Could not retrieve rigid body dof "<<j
                        <<" for population from index\n";
                ++nerr;
            } else {
                tmp=get_dof_info(uid);
                if (tmp.dof_kind!=rigid_body_xyz) {
                    prf::cerr<<"DOF kind mismatch with index for rigidbody "
                            <<"dof "<<j<<" of the population.\n";
                    ++nerr;
                }
                if (tmp.specific_global_index!=j) {
                    prf::cerr<<"Serial number mismatch with index for "
                            <<"rigidbody dof "<<j<<" of the population.\n";
                    ++nerr;
                }
            }
        }
        //backbone
        int nbdof=0;
        for (size_t i=0;i<chv.size();++i) nbdof+=Chain(i)->numBBdof();
        for (int j=0;j<nbdof;++j) {
            int uid=dofindex.get_uid_from_type_index(backbone_torsion_angle,j);
            if (uid<0) {
                prf::cerr<<"Could not retrieve backbone dof "<<j
                        <<" for population from index\n";
                ++nerr;
            } else {
                tmp=get_dof_info(uid);
                if (tmp.dof_kind!=backbone_torsion_angle) {
                    prf::cerr<<"DOF kind mismatch with index for backbone "
                            <<"dof "<<j<<" of the population.\n";
                    ++nerr;
                }
                if (tmp.specific_global_index!=j) {
                    prf::cerr<<"Serial number mismatch with index for "
                            <<"backbone dof "<<j<<" of the population.\n";
                    ++nerr;
                }
            }
        }
        //sidechain
        int nsdof=0;
        for (size_t i=0;i<chv.size();++i) nsdof+=Chain(i)->numRTdof();
        for (int j=0;j<nsdof;++j) {
            int uid=dofindex.get_uid_from_type_index(sidechain_torsion_angle,j);
            if (uid<0) {
                prf::cerr<<"Could not retrieve sidechain dof "<<j
                        <<" for population from index\n";
                ++nerr;
            } else {
                tmp=get_dof_info(uid);
                if (tmp.dof_kind!=sidechain_torsion_angle) {
                    prf::cerr<<"DOF kind mismatch with index for sidechain "
                            <<"dof "<<j<<" of the population.\n";
                    ++nerr;
                }
                if (tmp.specific_global_index!=j) {
                    prf::cerr<<"Serial number mismatch with index for "
                            <<"sidechain dof "<<j<<" of the population.\n";
                    ++nerr;
                }
            }
        }
        if (nerr!=0) {
            prf::cerr<<"Population> "<<nerr<<" errors in DOF index.\n";
            return 0;
        } else {
            Logger()(10)<<"Population> DOF index is clean.\n";
            return 1;
        }
    }
}
