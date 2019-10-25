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

#include <MC.hh>
#include <ProgUtils.hh>
#include <Bias.hh>
#include <HydrogenBond.hh>
#include <Hydrophobicity.hh>
#include <LocExVol.hh>
#include <ExVol.hh>
#include <Rot.hh>
#include <Pivot.hh>
#include <BGS.hh>
#include <RCUpdate.hh>
#include <Translation.hh>
#include <Rotation.hh>
#include <UpdateProbs.hh>
#include <RandomNumber.hh>
using namespace prf;
void help_function();
int main(int argc, char *argv[]) 
{
    //First list out the options you want to handle on command line
    ProgArgs prog;
    //Format: long name, short name, number of arguments
    prog.option("sequence","s",1); 
    prog.option("number_of_chains","nc",1);
    prog.option("init_type","i",1);
    prog.option("temperature","t",1);
    prog.option("seed","d",1);
    prog.option("number_of_cycles","ncyc",1);
    prog.option("length_of_cycles","lcyc",1);
    prog.option("write_interval","nrt",1);
    prog.option("pdb_write_interval","npdb",1);
    prog.option("update_probability_file","upfile",1);
    prog.option("log_level","ll",1);
    prog.option("help","h",0);
    //Analyze command line to see what options are passed
    prog.analyze(argc,argv);
    //Complain, if no arguments are passed
    if (argc==1) {
	prf::cout<<"Usage:\n\n"
		<<argv[0]<<" [OPTIONS] \n\n"
		<<"where, OPTIONS means one or more of the following\n"
		<<"in arbitrary order.\n";
	prog.write_available();
	return 0;
    }
    //Call the help function and quit if "-h" option is given
    if (prog.option_given("h")) {
	help_function();
	return 0;
    }
    //Define variables you want to modify and their default values
    string sequence="*GAVLISTCMPDNEQHKRFYW*";
    string init_type="random";
    double T=300*kelvin_in_pru;
    long seed=1111,ncyc=10000,nrt=1000,npdb=1000;
    int nc=1,lcyc=-1,loglevel=20;
    string update_probability_file="updateprobs.dat";
    //Modify the defaults according to the command line
    if (prog.option_given("s")) sequence=prog.option("s");
    if (prog.option_given("nc")) nc=atoi(prog.option("nc").c_str());
    if (prog.option_given("i")) init_type=prog.option("i");
    if (prog.option_given("t"))
	T=strtod(prog.option("t").c_str(),NULL)*kelvin_in_pru;
    if (prog.option_given("d")) seed=atol(prog.option("d").c_str());
    if (prog.option_given("ncyc")) ncyc=atol(prog.option("ncyc").c_str());
    if (prog.option_given("lcyc")) lcyc=atoi(prog.option("lcyc").c_str());
    if (prog.option_given("nrt")) nrt=atol(prog.option("nrt").c_str());
    if (prog.option_given("npdb")) npdb=atol(prog.option("npdb").c_str());
    if (prog.option_given("upfile")) 
	update_probability_file=prog.option("upfile");
    if (prog.option_given("ll")) loglevel=atoi(prog.option("ll").c_str());
    prf::Logger::verbosity=loglevel;
    //Declare known energies in PROFASI
    Bias bias;
    HBMM hbmm;
    HBMS hbms;
    Hydrophobicity hp;
    LocExVol lexv;
    ExVol exv;
    //Declare known updates in PROFASI
    Rot rot;
    Pivot pivot;
    BGS bgs;
    RCUpdate rcu;
    Translation trans;
    Rotation rotn;
    //Set up a random number generator
    RandomNumber ran;
    if (prog.option_given("d")) ran.setSeed(seed); else ran.RandomizeState();
    //Initialize library of groups, geometric parameters...
    Groups::initGroups();
    //Set up the population of proteins
    Population p;
    p.RandomNumberGenerator(&ran);
    p.AddProtein(sequence,nc);
    p.Init(); 
    p.InitCoord(init_type);
    p.EnforceBC();
    AtomCoordinates::update(0,p.NumberOfAtoms());
    //Why not save the initial conformation as a PDB !
    FILE *curpdb=fopen("initial.pdb","w");
    p.WritePDB(curpdb);
    fclose(curpdb);
    //Set up a Markov Chain Monte Carlo Simulator
    MC mc;
    mc.EnergyTerm(&bias);
    mc.EnergyTerm(&hbmm);
    mc.EnergyTerm(&hbms);
    mc.EnergyTerm(&hp);
    mc.EnergyTerm(&lexv);
    mc.EnergyTerm(&exv);
    mc.AddUpdate(&rot);
    mc.AddUpdate(&pivot);
    mc.AddUpdate(&bgs);
    mc.AddUpdate(&rcu);
    if (nc>1 || prog.option_given("frg")) {
	mc.AddUpdate(&trans);
	mc.AddUpdate(&rotn);
    }
    mc.Connect(&p);
    mc.temperature(T);
    mc.SetCycleLength(lcyc);
    mc.RandomNumberGenerator(&ran);
    mc.ReadProb(update_probability_file);
    mc.Setup();
    mc.InitEnergies();
    MCTotEnergy etot;
    etot.assoc(&mc);
    etot.refresh();
    //Now the main loop, the MCMC evolution
    //Open a run-time history file
    Output rt("rt");
    string sep="  ";
    double mine=etot.Value();
    for (long icyc=0;icyc<ncyc;++icyc) {
	if (icyc%nrt==0) 
	    prf::cout<<"cycle number "<<icyc<<" **************\n";
	mc.RunCycle();
	mc.EReset();
	etot.refresh();
	if (icyc%npdb==0) {
	    curpdb=fopen("current.pdb","w");
	    fprintf(curpdb,"REMARK   CYCLE=%li   ENERGY=%lf\n",icyc,etot.Value());
	    p.WritePDB(curpdb);
	    fclose(curpdb);
	    if (etot.Value()<mine) {
		p.WritePDB("minen.pdb");
		mine=etot.Value();
	    }
	}
	if (icyc%nrt==0) {
	    rt<<icyc<<sep;
	    rt<<etot.Value()<<sep;
	    bias.refresh();
	    hbmm.refresh();
	    hbms.refresh();
	    hp.refresh();
	    lexv.refresh();
	    exv.refresh();
	    rt<<bias.Value()<<sep<<hbmm.Value()<<sep<<hbms.Value()<<sep
		    <<hp.Value()<<sep<<lexv.Value()<<sep<<exv.Value()<<"\n";
	}
    }
    rt.close();
    mc.output_statistics("updates.stat");
    return 0;
}
//Help function for more detailed description of the options and 
//their meanings...
void help_function() 
{
    
}
