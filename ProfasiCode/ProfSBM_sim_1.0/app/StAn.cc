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

#include "StAn.hh"
#include <Aux/PDBReader.hh>
#include <Aux/RMSD_Utils.hh>
#include <list>
#include <Aux/ProgUtils.hh>

using std::list;
using namespace prf;
using namespace prf_utils;

StAn::StAn()
{
    pre_init();
}

StAn::~StAn()
{

}

void StAn::pre_init()
{
    Groups::initGroups();
    AminoAcid::initCommon();
    H.unsetSwitch("histograms");
    H.setPopulation(&p);
    H.set_block_props(1,"temperature");
    AtomCoordinates::SetBox(500);
}

void StAn::re_init_obs()
{
    H.re_init_obs();
}

bool StAn::composition_changed()
{
    bool changed=false;
    if (p.NumberOfChains()!=thecomposition.size()) changed=true;
    else {
        for (size_t ic=0;(!changed)&&(ic<thecomposition.size());++ic) {
            if (thecomposition[ic].size()!=p.Chain(ic)->numLigands()) changed=true;
            else {
                for (size_t jr=0;(!changed)&&(jr<thecomposition[ic].size());++jr) {
                    changed=(thecomposition[ic][jr]!=p.Chain(ic)->memberLigand(jr)->OLC());
                }
            }
        }
    }
    return changed;
}

int StAn::load(std::string ifile)
{
    Logger blog;
    size_t icolon=ifile.find(':');
    std::string sel, filnm;
    if (icolon<ifile.size()-1) sel=std::string(ifile,icolon+1);
    else sel="1:*";
    filnm=std::string(ifile,0,icolon);

    std::vector<std::string> iparts;
    split_str(filnm,'.',iparts);
    std::string iformat=iparts.back();

    if (iformat=="xml") {
        prf_xml::XML_Node *root=prf_xml::get_xml_tree(ifile);
        if (root!=NULL and root->child("population")!=NULL)
            p.Read_XML(root->child("population"));
        else {
            prf::cerr<<"No population found in "<<ifile<<"\n";
            return 1;
        }
    } else if (iformat=="pdb") {
        PDBReader pdb;
        pdb.set_file(filnm);

        if (pdb.read_matrix()==0) {
            prf::cerr<<"Failed to read in data from "<<filnm<<"\n";
            return 1;
        }

        pdb.check_file();
        blog(50)<<"input file : "<<filnm<<" selection "<<sel<<"\n";
        list<SelRes> selection;
        pdb.mk_selection(sel,selection);
        blog(50)<<"Initial selection size = "<<selection.size()<<"\n";

        p.AddProtein(selection);
        list<AtomRecord> rcd;
        pdb.records(selection,rcd);

        p.Init();
        p.InitCoord("stretched");
        std::vector<bool> specified(p.NumberOfAtoms(),false);
        p.ImportStructure(rcd,specified);
        p.guess_missing_coordinates(specified);
    } else {
        prf::cerr<<"Unknown input format.\n";
        return 1;
    }
    p.Reconstruct();
    p.EnforceBC();
    if (composition_changed()) re_init_obs();

    return 0;
}

void StAn::save(std::string ofile)
{
    size_t extn=ofile.find_last_of('.');
    std::string oformat=ofile.substr(extn+1);
    int infmt=0, curTindex=0;
    unsigned long icyc=0;
    if (oformat=="pdb") infmt=1;
    else if (oformat=="xml") infmt=2;
    else if (oformat=="tcn" or oformat=="tconf" or oformat=="tcnf") infmt=3;
    else if (oformat=="bcn") infmt=4;
    else {
        prf::cerr<<"Can not determine output format from file name "<<ofile<<"\n";
        return;
    }
    double etot=0;
    if (H.get_obs("Etot")!=NULL) etot=H.get_obs("Etot")->Value();
    p.SaveSnapshot(infmt,ofile,icyc, curTindex,etot);
}

Matrix<double> StAn::profile(std::string exprs,std::string rng)
{
    Matrix<double> ans;
    return ans;
}

void StAn::interactive_session()
{
    InstructionString s;
    std::string line;

    do {
        if (s.type()!=INCOMPLETE) prf::cout<<"stan> ";
        else prf::cout<<"continuation> ";
        getline(std::cin,line);
        if (s.type()==INCOMPLETE) s.append(line); else s.str(line);
        if (s.type()==COMMAND) {
            process_command(s);
        }
    } while (s.type()!=QUIT);
}

void StAn::process_command(InstructionString s)
{
    Observable *o=NULL;
    if (s.head()=="load") {
        if (s.n_parts()==2) {
            load(s.part(1));
        }
    } else if (s.head()=="save") {
        if (s.n_parts()==2) {
            save(s.part(1));
        }
    } else if (s.head()=="list_obs") {
        for (int i=0;i<H.num_obs();++i) {
            prf::cout<<i<<"\t"<<H.get_obs(i)->Name()<<"\n";
        }
    } else if (s.head()=="log_level") {
        prf::Logger::verbosity=atoi(s.tail().str().c_str());
    } else if (s.head()=="show" or s.head()=="print") {
        show(s.tail());
    } else if ((o=H.get_obs(s.head()))!=NULL) {
        o->refresh();
        prf::cout<<o->Value()<<"\n";
    } else if (!H.parseCommand(s)) {
        prf::cout<<"Unknown command "<<s.head()<<"\n";
    }
}

void StAn::show(InstructionString s)
{
    Observable *o=H.get_obs(s.head());
    if (o!=NULL) {
        o->refresh(0);
        prf::cout<<o->Name()<<" = "<<o->Value()<<"\n";
    } else prf::cout<<s.head()<<": Unknown observable.\n";
}

void StAn::parseCommands(std::string flnm)
{
    if (TestFile_r(flnm)) {
        H.ParseCommands(flnm,".");
    }
}

int main(int argc, char *argv[])
{
    prf_utils::ProgArgs optn;
    optn.option("settings_file","s",1,"(Set up observables from a given settings file)");
    optn.analyze(argc,argv);
    StAn interface;
    if (optn.n_spare_args()==1) interface.load(optn.spare_args(0));
    if (optn.option_given("s")) {
        interface.parseCommands(optn.option("s"));
    }
    interface.interactive_session();
    return 0;
}
