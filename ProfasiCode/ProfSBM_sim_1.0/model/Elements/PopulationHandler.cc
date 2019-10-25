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

#include "PopulationHandler.hh"
#include "../Aux/prf_xml.hh"
#include "../Aux/PDBReader.hh"
#include <list>

PopulationHandler::PopulationHandler() : HandlerBase()
{
    initstate="random";
    usrspc_box_size=false;
    AminoAcid::initCommon();
    par.option("box_length","box_length",1,"(Size of periodic box)");
    par.option("add_chain","ac",2,"(number of chains and sequence in quotes)");
    par.option("add_chain_pdb","acpdb",2,"(number of chains and a pdb file with a selection)");
    par.option("set_population","setpop",1,"(name of the xml file with the population info)");
    par.option("init_config","init_config",1,"(\"stretched\" or \"random\" or \"file:://filename\")");
    par.option("cis","cis",2,"(which chain and which peptide bond)");
    par.option("void_ends","vd",0,"(Leave chain ends uncharged in absense of end groups)");
}
/**
\page settings_population Commands to set up the population
\li \b --add_chain or \b -ac : Add chains to the population. Example: --add_chain 3
"<Acetyl*KLVFFAE*NH2>" <br> The sequence has to be given in the format explained in
\ref seq_input. The angular brackets in this example can be dropped. They were
necessary for ProFASi versions before 1.5. Also, note that this command can be
repeated on the command line or the settings file to add more chains. Something
like <br> -ac 2 "*KFFE AAAK KFFE*" -ac 2 "*KKAFAFAFEE*" <br>
would add 4 chains.
\li \b --add_chain_pdb or \b -acpdb : Add a sequence from a PDB file. Example:
--add_chain_pdb 1 1GB1.pdb::A,41,56<br> The sequence can be a selection from the
PDB file, specified by the ProFASi selection rules (See \ref prf_sel_fils).
\li \b --set_population or \b -setpop : Import sequence and structure of all
chains from a ProFASi XML structure file. Example: -setpop minen.xml
\li \b --init_config or \b -init_config : Initial configuration. Possible
 are "random", "stretched" or "file://a_text_config_file"
 \li \b --cis or \b -cis : Declare a certain peptide bond to be "cis". The command
 takes two arguments, the chain id and the peptide bond id. Example: -cis 0 15
 \li \b --void_ends or \b -vd : Even if there are no capping groups, leave the
 chain ends uncharged
 \li \b --box_length or \b -b : Size of the periodic box, in Angstroems
*/
PopulationHandler::~PopulationHandler() {}

int PopulationHandler::parseCommand(InstructionString s)
{
    Logger(30)<<"PopulationHandler> processing command "<<s.str()<<"\n";
    if (s.head()=="box_length") boxSize(strtod(s.part(1).c_str(),NULL));
    else if (s.head()=="add_chain") {
        if (s.n_parts()<3) {
            prf::cerr<<"Incorrect syntax for add_chain. Discarded. Example usage:\n";
            prf::cerr<<"add_chain 1 < Acetyl * KLVFFAE * NH2 >\n";
        } else {
            p.AddProtein(s.tail().tail().str(),atoi(s.part(1).c_str()));
        }
    } else if (s.head()=="add_chain_pdb") {
        if (s.n_parts()<3) {
            prf::cerr<<"Incorrect syntax for add_chain_pdb. Discarded.\n";
        } else {
            p.AddProtein(atoi(s.part(1).c_str()),s.part(2));
        }
    } else if (s.head()=="set_population") {
        if (s.n_parts()<2) {
            prf::cerr<<"Incorrect syntax for set_population. Discarded.\n";
        } else {
            std::string fullnm=s.part(1);
            size_t idot=(fullnm.size()),icolon=(fullnm.size());
            idot=fullnm.find_last_of('.');
            if (idot<fullnm.size()) icolon=fullnm.find(':',idot);
            std::string extn=fullnm.substr(idot+1,(icolon-idot-1));
            std::string flnm=fullnm.substr(0,icolon);
            std::string sel=fullnm.substr(icolon+1);
            if (extn=="xml") read_xml_pop(s.part(1));
            else if (extn=="pdb") read_pdb_pop(flnm,sel);
            else {
                prf::cerr<<"set_population: Import is only available for file "
                        <<"extensions xml and pdb. Parsing \""<<s.str()
                        <<"\", the full file name was inferred to be \""
                        <<flnm<<"\", with extension \""<<extn<<"\"";
                if (!sel.empty()) prf::cerr<<", and selection "<<sel;
                prf::cerr<<"\n";
            }
        }
    } else if (s.head()=="init_config") {
        if (s.n_parts()<2) {
            prf::cerr<<"Incorrect syntax for init_config. Discarded.\n";
        } else {
            initstate=s.tail().str();
        }
    } else if (s.head()=="cis") {
        p.setCis(atoi(s.part(1).c_str()),atoi(s.part(2).c_str()));
    } else if (par.option_given("vd")) p.charged_ends(false,false);

    return 1;
}

void PopulationHandler::setupPeriodicBox()
{
    double ans;
    Logger blog;

    if (usrspc_box_size) ans=boxl;
    else {
        ans=3.8*(p.LongestChain()->numLigands())+4.5;
        ans*=pow(p.NumberOfChains(),1.0/3);
        blog(3)<<"PopulationHander> Using automatically calculated box size "
                <<ans<<" Angstroms\n";
    }

    AtomCoordinates::SetBox(ans);

}

void PopulationHandler::boxSize(double xx, bool usbit)
{
    boxl=xx;usrspc_box_size=usbit;
}

int PopulationHandler::init_pop()
{
    p.Init();
    if (not p.initialized()) return 0;
    setupPeriodicBox();
    return 1;
}

int PopulationHandler::init_coords()
{
    int initres=p.InitCoord(initstate);
    if (initstate!="random" && initstate!="stretched" && initres==1) return 1;
    else return 0;
}

int PopulationHandler::reconstruct()
{
    p.EnforceBC();
    p.Reconstruct();
    AtomCoordinates::update(0,p.NumberOfAtoms());
    return 1;
}

int PopulationHandler::read_xml_pop(std::string xmlfile)
{
    Logger blog;
    prf_xml::XML_Node * root = prf_xml::get_xml_tree(xmlfile);
    if (root==NULL) {
        prf::cerr<<"PopulationHandler> Could not retrieve a valid top level XML "
            <<"node from "<<xmlfile<<"\n";
        return 0;
    } else {
        blog(15)<<"PopulationHandler> Retrieved XML node with "
            <<root->n_children()<<" sub-nodes\n";
        for (size_t i=0;i<root->n_children();++i) {
            blog<<root->child(i)->name()<<"\n";
        }
    }
    prf_xml::XML_Node * popl=NULL;
    if (root->name()=="population") popl=root;
    else popl=root->child("population");
    if (popl==NULL) {
        prf::cerr<<"PopulationHandler> No population found in "<<root->name()<<"\n";
        return 0;
    }
    p.Read_XML(popl);
    if (root->child("box_length")!=NULL) {
        boxSize(strtod(root->child("box_length")->value().c_str(),NULL));
    }
    setupPeriodicBox();
    initstate="none"; //don't assign anything later
    reconstruct();
    delete root;
    return 1;
}

int PopulationHandler::read_pdb_pop(std::string flnm,std::string sel)
{
    PDBReader pdb(flnm.c_str());

    if (pdb.read_matrix()==0) {
        prf::cerr<<"Could not read in information from "<<flnm<<"\n";
        return 0;
    }

    std::list<SelRes> slct;
    pdb.mk_selection(sel,slct);
    p.clear();
    p.AddProtein(slct,1);
    if (!init_pop()) return 0;
    std::list<AtomRecord> rcds;
    pdb.records(slct,rcds);
    std::vector<bool> crdchk(p.NumberOfAtoms(),false);
    p.ImportStructure(rcds,crdchk,0);
    int unspec=std::count(crdchk.begin(),crdchk.end(),false);
    if (unspec!=0)
        prf::cerr<<"PopulationHandler> Number of unspecified "
                <<"coordinates in file "<<flnm<<" = "<<unspec<<"\n";
    p.guess_missing_coordinates(crdchk);
    for (int ich=0;ich<p.num_chains();++ich)
        p.Chain(ich)->calc_torsions(crdchk);
    initstate="none";
    return 1;
}
