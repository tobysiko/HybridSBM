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
#include <Aux/fileutils.hh>
#include <Aux/PDBReader.hh>
#include <Elements/Population.hh>
#include <Observables/ContactOrder.hh>

using namespace prf;

using namespace prf_utils;
using std::string;
using std::vector;
using std::list;

int main(int argc, char *argv[])
{
    ProgArgs par;
    par.option("log_level","ll",1);
    par.option("cutoff","c",1);
    par.new_switch("verbose","v",false);
    par.analyze(argc,argv);

    if (argc==1) {
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [parameters] pdbfile:model:selections\n\n";
        prf::cout<<"parameters could be any number and combination of ...\n";
        par.write_available();
        prf::cout<<"Example: \n\n";
        prf::cout<<argv[0]<<" -l 3 -c 5.0 minen.pdb::A\n\n";
        prf::cout<<"Contact order is evaluated after initializing a population\n"
                 <<"based on the sequence and structure taken from the selection \n"
                 <<"provided in the command line. It only makes sense for contiguous\n"
                 <<"selections on the same chain. So, using the full flexibility of\n"
                 <<"PROFASI pdb selection rules with this program will most likely \n"
                 <<"lead to meaningless answers.\n\n";
        return 0;
    }

    int loglev=10;
    prf::clog.open(".profasi_logfile","w");

    if (par.option_given("ll")) loglev=atoi(par.option("ll").c_str());

    prf::Logger::verbosity=loglev;
    double dcut=6.0;

    if (par.option_given("c")) dcut=strtod(par.option("c").c_str(),NULL);

    bool verbose=(par.switch_given("v"));

    string sel, filnm,ifile="inputfile.pdb";

    if (par.n_spare_args() >=1) ifile=string(par.spare_args(0));

    size_t icolon=ifile.find(':');

    if (icolon<ifile.size()-1) sel=string(ifile,icolon+1);
    else sel="1:*";

    filnm=string(ifile,0,icolon);
    Groups::initGroups();
    AminoAcid::initCommon();
    PDBReader pdb;
    pdb.set_file(filnm);

    if (pdb.read_matrix()==0) {
        prf::cerr<<"Failed to read in data from "<<filnm<<"\n";
        return 1;
    }

    if (verbose) prf::cout<<"input file : "<<filnm<<" selection "<<sel<<"\n";

    list<SelRes> selection;
    list<AtomRecord> rcd;
    pdb.mk_selection(sel,selection);
    Population p;
    p.AddProtein(selection);
    p.Init();
    p.InitCoord("stretched");
    pdb.records(selection,rcd);

    if (verbose) prf::cout<<"Selection has total size "<<selection.size()
                              <<" residues\n";

    if (verbose) prf::cout<<"That contains "<<rcd.size()<<" atom records\n";

    std::vector<bool> specified(p.NumberOfAtoms(),false);
    p.ImportStructure(rcd,specified);
    p.guess_missing_coordinates(specified);
    FILE *fp=fopen("new.pdb","w");
    p.WritePDB(fp);
    fclose(fp);

    ContactOrder CO;
    CO.set_cutoff(dcut);
    CO.setPopulation(&p);
    CO.init_obs();
    CO.refresh();

    if (verbose) prf::cout<<"Relative contact order = ";

    prf::cout<<CO.Value()<<"\n";
    return 0;
}

