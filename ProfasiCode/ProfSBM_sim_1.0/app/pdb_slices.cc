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

#include <Aux/ProgUtils.hh>
#include <Aux/profasi_io.hh>
#include <Aux/PDBReader.hh>
#include <Aux/fileutils.hh>
#include <Aux/RMSD_Utils.hh>
#include <Elements/GroupLib.hh>

using namespace prf;

using namespace prf_utils;
using std::string;
using std::list;
using std::vector;

int main(int argc, char *argv[])
{
    ProgArgs par;
    par.option("logfile","l",1);
    par.option("verbose_mode","v",0);
    par.option("using","u",1);
    par.option("output_file","o",1);
    par.option("reconstruct","r",1);
    par.option("log_level","ll",1);
    par.analyze(argc,argv);

    if (argc==1) {
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] file.pdb\n\n";
        prf::cout<<"Options could be any number and combination of ...\n";
        par.write_available();
        prf::cout<<"To extract the 15th model from a PDB file containing \n"
        <<"many models, write\n\n";
        prf::cout<<argv[0]<<" file.pdb:15 -o file_15.pdb\n\n";
        prf::cout<<"To extract the chain B of model 5, write\n\n";
        prf::cout<<argv[0]<<" file.pdb:5:B -o file_5B.pdb\n\n";
        prf::cout<<"To extract the residues 21 to 52 of chain B\n"
        <<"of model 5, write\n\n";
        prf::cout<<argv[0]<<" file.pdb:5:B,21,52 -o file_5B_21_52.pdb\n\n";
        prf::cout <<"To make a pdb file out of the backbone and CB atoms of \n"
        <<"a file, write\n\n"
        <<argv[0]<<" file.pdb -u \"+BB+CB\" -o file2.pdb\n\n"
        <<"To renumber atom and residue indices in the output pdb file, use \n"
        <<argv[0]<<" file.pdb -u \"+BB+CB\" -o file2.pdb -r \"atom_index&residue_index\"\n\n";
        prf::cout<<"To supress or increase the amount of log messages written\n"
        <<"in the log file, use the -ll option to set a log level.\n"
        <<"Higher values mean more messages. Default is 10.\n\n";
        return 0;
    }

    // Some initialisation
    Logger blog;

    string lgfl=".profasi_logfile",opfl="pdb_slices.output";

    int loglvl=10;

    if (par.option_given("l")) lgfl=par.option("l");

    if (par.option_given("ll")) loglvl=atoi(par.option("ll").c_str());

    if (par.option_given("o")) opfl=par.option("o");

//    prf::clog.open(lgfl.c_str());

    bool verbose_mode=false;

    if (par.option_given("v")) verbose_mode=true;

    prf::Logger::verbosity=loglvl;

    if (par.n_spare_args()<1) {
        prf::cerr<<"Needs a pdb file to act upon.\n";
        return 1;
    }

    prf::Groups::initGroups();

    string filename,myfilter;
    filename=par.spare_args(0);
    size_t icolon=filename.find(':');
    string filenm(filename,0,icolon);
    string selections;

    if (icolon<filename.size()-1) selections=string(filename,icolon+1);
    else selections="1:*";

    if (verbose_mode) prf::cout<<"File name :"<<filenm
        <<", selections :"<<selections<<"\n";

    PDBReader pdb(filenm.c_str());

    if (pdb.read_matrix()==0) {
        prf::cerr<<"Could not read in information from "<<filenm<<"\n";
        return 0;
    }

    blog(8)<<"Number of Chains = "<<pdb.num_chains()<<"\n";

    for (int i=0;i<pdb.num_chains();++i) {
        blog(8)<<"Chain "<<i<<"="<<pdb.chain_name(i)<<"\n";
    }

    list<SelRes> slct;

    pdb.mk_selection(selections,slct);
    blog(8)<<"Selection has "<<slct.size()<<" entries\n";
    list<AtomDescriptor> lst;
    pdb.descriptors(slct,lst);

    myfilter="BB";

    if (par.option_given("u")) myfilter=par.option("u");

    RMSD_Utils utils;

    utils.make_filter(myfilter);

    utils.apply_filter(lst);

    vector<int> slctd(lst.size(),0);

    int j=0;

    for (list<AtomDescriptor>::iterator i=lst.begin();i!=lst.end();++i) {
        slctd[j++]=i->int_label;
    }

    list<AtomRecord> rcd;

    pdb.export_records(slctd,rcd);

    if (par.option_given("r")) {
        if (verbose_mode) prf::cout<<"The PDB info will be regenerated.\n";

        string rec_str=par.option("r");

        vector<string> recparts;

        split_str<vector<string> >(rec_str,'&',recparts);

        bool recatno=false, recresno=false;

        int recatst=1,recrsst=atoi(rcd.begin()->descriptor().ires.c_str());

        for (size_t ri=0;ri<recparts.size();++ri) {
            if (recparts[ri].find(string("atom_index"))<recparts[ri].size()) {
                vector<string> atnostr;
                split_str<vector<string> >(recparts[ri],':',atnostr);

                if (atnostr.size()>1) recatst=atoi(atnostr[1].c_str());

                recatno=true;
            }

            if (recparts[ri].find(string("residue_index"))<recparts[ri].size()) {
                vector<string> rsnostr;
                split_str<vector<string> >(recparts[ri],':',rsnostr);

                if (rsnostr.size()>1) recrsst=atoi(rsnostr[1].c_str());

                recresno=true;
            }
        }

        if (recatno) {
            if (verbose_mode) prf::cout<<"Renumbering atoms.\n";

            int curat=recatst;

            for (list<AtomRecord>::iterator it=rcd.begin();it!=rcd.end();++it) {
                it->descriptor().iatom=curat++;
            }
        }

        if (recresno) {
            bool all_res_index_are_numbers=true;

            for (list<AtomRecord>::iterator it=rcd.begin();it!=rcd.end();++it) {
                string resin=(it->descriptor().ires);

                for (size_t irc=0;irc<resin.size();++irc) {
                    if (!isdigit(resin[irc])&&resin[irc]!=' ' &&resin[irc]!='-') {
                        all_res_index_are_numbers=false;
                        break;
                    }
                }

                if (!all_res_index_are_numbers) break;
            }

            if (all_res_index_are_numbers) {
                int iresstart=recrsst;

                if (verbose_mode) prf::cout<<"Renumbering residues offset by "
                    <<iresstart<<"\n";

                for (list<AtomRecord>::iterator it=rcd.begin();it!=rcd.end();++it) {
                    int newires=atoi(it->descriptor().ires.c_str())-iresstart+1;
                    char irres[4];
                    sprintf(irres,"%3d",newires);
                    it->descriptor().ires=string(irres);
                }
            } else {
                prf::cerr<<"Not all residue numbers in this file are positive "
                <<"integers.\n"
                <<"I was unable to reassign the residue numbering.\n";
            }
        }
    }

    for (list<AtomRecord>::iterator it=rcd.begin();it!=rcd.end();++it) {
        it->build_pdb_line_from_fields();
    }

    Output op;

    if (par.option_given("o")) {
        op.open(par.option("o").c_str());
    } else {
        op.Attach(stdout);
    }

    for (list<AtomRecord>::iterator it=rcd.begin();it!=rcd.end();++it) {
        op<<it->pdb_line()<<"\n";
    }

    return 0;
}

/**
\page pdb_slices Extracting parts of a PDB file
This is a utility to extract a portion of a given PDB file, like a particular model, one of the many chains, a certain range of residues in some chain, or only a certain kind of atoms. It works with the selection and filter rules as described in \ref prf_sel_fils . \n\n

Usage: \n\n
pdb_slices [OPTIONS] input.pdb:model:selections \n\n
\section Examples
<tt>pdb_slices abcd.pdb:5 -o model5.pdb</tt> \n\n
This would create a new pdb file model5.pdb with the 5th model in abcd.pdb.\n

<tt>pdb_slices abcd.pdb:5:A -o model5_chainA.pdb </tt>\n\n
This creates a new pdb file with the chain A of model 5 in abcd.pdb.\n

<tt>pdb_slices abcd.pdb::A,2,30 -u "+BB+CB" -o bbcb.pdb</tt> \n\n
This uses the residues labeled 2 through 30 in chain A of the first model in abcd.pdb. Only backbone and C-beta atoms are present in the new pdb file. \n

<tt>pdb_slices abcd.pdb::A,2,30 -u "+BB" -r atom_index -o new.pdb</tt> \n\n
In this case, the atom numbers in the new pdb file will be 1,2,3..., instead of whatever numbers the selected atoms had in the original file. Similarly, to have residues renumbered from 1, replace atom_index with residue_index. To have both renumbered, use "-r atom_index\&residue_index".

*/
