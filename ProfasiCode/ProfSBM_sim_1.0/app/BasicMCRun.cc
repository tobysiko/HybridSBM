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

#ifdef PARALLEL
#include <mpi.h>
#endif
#include "BasicMCRun.hh"
#include <Aux/Timer.hh>
#include <Aux/profasi_version.hh>

using namespace prf;
using namespace UnivConstants;
using std::min;
using std::string;

BasicMCRun::BasicMCRun()
{
    progheader="Canonical Metropolis";
    swtch["stdout_redirect"]=false;
    swtch["stderr_redirect"]=false;
    swtch["preliminary_relaxation"]=true;
    swtch["thermalisation"]=true;
    swtch["limited_time"]=false;
    swtch["resume_mode"]=false;
    // resume_mode is true when starting from the middle of a previous run
    swtch["continue"]=true;
    // continue is true if it starts from the end of the last run. Also,
    // if there is no previous run, it is true.

    for (map<string,bool>::iterator it=swtch.begin();
    it!=swtch.end();++it) {
        usrspc[it->first]=false;
    }
    
    debug = false;
    
    log_level=10;
    prf::Logger::verbosity=log_level;
    NRT=1000;

    NTHERM=25000;
    ICONF=1000;
    IAVG=10000;
    MCCYC=100000000;
    RMCCYC=(unsigned int) (-1);
    INORM=1;
    NTESTI=10000;
    rem_progress_writes=10;
    ncyc_per_T_updt=1;
    locdir=".";
    pdbfil="current";minpdb="minen";
    statfile="updates.stats";
    progname="BasicMCRun";
    Tcur=0.46;
    Trelax=1.0;
    curTindex=0;
    NTMP=1;
    icyc=icyc0=0;
    nblk=0;
    mc=&canon;
    available_run_time=1e20;
    myrank=0;nruns=1;
    snapshot_format=1;
    
    
    optn.option("num_cycles","ncyc",1,"(number of MC sweeps (Note:sweeps!=steps))");
    optn.option("num_therm_cycles","ntherm",1,
                "(number of thermalisation sweeps)");
    optn.option("available_time","time",1,
                "(time allocated for the run)");
    optn.option("rt_write_freq","nrt",1,
                "(interval for writing to the run-time history)");
    optn.option("conf_write_freq","iconf",1,
                "(interval of saving state to conf file)");
    optn.option("num_progress_reports","nreports",1,
                "(number of reports about projected finishing time etc. )");
#ifndef PARALLEL
    optn.option("rank","r",1,
                "(output directory number, if not controlled by MPI)");
    optn.option("nruns","nr",1,
                "(if running multiple runs, but without MPI)");
#endif
    optn.option("log_level","ll",1,"(verbosity level of log messages)");
    optn.new_switch("continue","continue",true,
                    "(Resume run from the end of a previous run (default))");
    optn.new_switch("use_relaxation_cycles","relax",true,
                "(relaxation sweeps before actual simulation (default))");
    optn.option("T_update_interval","ntup",1,
                "(Interval at which temperature update is attempted)");
    optn.new_switch("help","h",false,"(ask for usage information)");
    optn.new_switch("debug","dbg",false,"(show debugging info)");
    optn.disable("T_update_interval");
    optn.settings_use_type(1);
}

int BasicMCRun::parseCommands(std::list<InstructionString> &cmds)
{
    for (std::list<InstructionString>::iterator it=cmds.begin();
    it!=cmds.end();++it) {
        parseCommand(*it);
    }
    return 1;
}

int BasicMCRun::parseCommand(InstructionString s)
{	
	if (s.head()=="dbg" || s.head()=="debug"){
		if (s.part(1)=="on"){
			debug = true;
			swtch["debug"]=true;
		}
		else debug = false;
		
		std::cout<<"debug?"<<s.part(1)<<" --> "<<debug<<"\n";
	}else if (s.head()=="nrt" || s.head()=="rt_write_freq")
        NRT=strtoul(s.part(1).c_str(),NULL,10);
    else if (s.head()=="num_therm_cycles" || s.head()=="ntherm")
        NTHERM=strtoul(s.part(1).c_str(),NULL,10);
    else if (s.head()=="use_relaxation_cycles") {
        swtch["preliminary_relaxation"]=(s.part(1)=="on"||s.part(1)=="true");
    } else if (s.head()=="conf_write_freq" || s.head()=="iconf")
        ICONF=strtoul(s.part(1).c_str(),NULL,10);
    else if (s.head()=="iavg" || s.head()=="avg_write_freq")
        IAVG=strtoul(s.part(1).c_str(),NULL,10);
    else if (s.head()=="ncycles" || s.head()=="nsweeps" ||
             s.head()=="mccyc" || s.head()=="num_cycles")
        MCCYC=strtoul(s.part(1).c_str(),NULL,10);
    else if (s.head()=="snapshot_format")
        snapshot_format=atoi(s.part(1).c_str());
    else if (s.head()=="T_update_interval") {
        ncyc_per_T_updt=atoi(s.tail().str().c_str());
        if (ncyc_per_T_updt<=0) ncyc_per_T_updt=1;
    } else if (s.head()=="log_level") {
        log_level=atoi(s.part(1).c_str());
        prf::Logger::verbosity=log_level;
        prf::cout<<progname<<"> Set logger verbosity level to "
                <<log_level<<"\n";
    } else if (s.head()=="available_time") {
        available_run_time=make_seconds(s.tail().str());
        swtch["limited_time"]=true;
    } else if (s.head()=="resume") {
        trajfile=s.part(1);
        RMCCYC=strtoul(s.part(2).c_str(),NULL,10);
        trajfile=locdir+"/"+trajfile;
        swtch["resume_mode"]=true;
    } else if (s.head()=="num_progress_reports")
        rem_progress_writes=(unsigned) atoi(s.tail().str().c_str());

    for (map<string,bool>::iterator sw=swtch.begin();sw!=swtch.end();++sw) {
        if (s.head()==sw->first) {
            if (s.part(1)=="on"||s.part(1)=="true") {
                swtch[s.head()]=true;usrspc[s.head()]=true;break;
            } else if (s.part(1)=="off"|| s.part(1)=="false") {
                swtch[s.head()]=false;usrspc[s.head()]=true;break;
            }
        }
    }

    return 1;
}

void BasicMCRun::disable(std::string anoption)
{
    optn.disable(anoption);
}

void BasicMCRun::show_help()
{
    optn.write_available();
}

BasicMCRun::~BasicMCRun() {}

int BasicMCRun::update_T()
{
    return 0;
}

bool BasicMCRun::should_suspend(bool myrec)
{
    return myrec;
}

void BasicMCRun::set_index(int i, unsigned long j) {}

void BasicMCRun::write_MC_averages()
{
    mc->updates_hander()->output_statistics(statfile);
}

void BasicMCRun::init_streams()
{
    std::string wrtmd=(swtch["continue"]?"a":"w");

    if (swtch["stdout_redirect"])
        prf::cout.open((mydir+"/stdout").c_str(), wrtmd.c_str());

    if (swtch["stderr_redirect"])
        prf::cerr.open((mydir+"/stderr").c_str(), wrtmd.c_str());

}

void BasicMCRun::writeTemperatures()
{
    Output fp((mydir+string("/temperature.info")).c_str());
    fp<<"# Index T(model)\t T(Kelvin)\tBeta\n";
    fp<<"0\t"<<Tcur<<"\t"<<Tcur*pru_in_kelvin<<"\t"<<(1.0/Tcur)<<"\n";
    fp.close();
}

int BasicMCRun::init_MC()
{
	mc->setDebug(debug);
    return mc->Setup();
}

void BasicMCRun::init_filenames()
{
    string ext="";
    switch (snapshot_format) {
        case 1: ext=".pdb"; break;
        case 2: ext=".xml"; break;
        case 3: ext=".tcn"; break;
        case 4: ext=".bcn"; break;
        default: ;
    };

    pdbfil=mydir+"/current"+ext;
    minpdb=mydir+"/minen"+ext;
    statfile=mydir+"/"+statfile;
    speedfile=mydir+"/execution_speed.dat";
    rnboosterfile=mydir+"/random_number_state";
    trajfile=mydir+"/traj";
    mc->set_output_prefix(mydir+"/");
}

void BasicMCRun::run_relaxation_cycles()
{
	
    int ncycl=10*(PH.population()->NumberOfChains());
    
    std::cout<< "Relaxing the system for "<<ncycl<<" Monte Carlo cycles "
    <<"to get rid of obvious steric clashes. \nRelaxation "
    <<"finishes before entering the main loop, and hence does "
    <<"not enter any measurements.\n";
    
    prf::clog<<"Relaxing the system for "<<ncycl<<" Monte Carlo cycles "
    <<"to get rid of obvious steric clashes. \nRelaxation "
    <<"finishes before entering the main loop, and hence does "
    <<"not enter any measurements.\n";
    mc->temperature(Trelax);

    for (int icy=0;icy<10*(PH.population()->NumberOfChains());++icy) {
        mc->RunCycle();
        ffh.interaction_potential()->reset_total();
    }

    mc->temperature(Tcur);



}


int BasicMCRun::Run()
{
    printMCconfig();
    H.printConfig();
    highlight();
    prf::cout<<"\n\n\n";

    if (swtch["preliminary_relaxation"]) run_relaxation_cycles();

    prf::cout<<"\nEnergies before going into main loop...\n";
    ffh.interaction_potential()->print_contributions(prf::cout);

    //The main loop ...


    PH.population()->SaveSnapshot(2, mydir+"/start.xml",icyc, curTindex,
                                  ffh.interaction_potential()->value());

    NTESTI=icyc+min((unsigned long) NTESTI,(MCCYC-icyc)/10);

    prf::cout<<"\n\nThis run will perform Monte Carlo cycles starting "
            <<"at count "<<icyc<<" up to count "<<MCCYC
            <<". This corresponds to "
            <<(MCCYC-icyc)<<" Monte Carlo cycles, or "
            <<((MCCYC-icyc)*mc->CycleLength())
            <<" elementary Monte Carlo updates.\n\n\n";

    if (swtch["thermalisation"] and icyc<NTHERM) {
        H.disableStatistics();
    }

    prf::clog.flush();prf::cerr.flush();prf::cout.flush();

    double emin(1e50);
    Timer timer;

    timer.space_axis_name("cycle");
    timer.set_new_goal(MCCYC);
    timer.n_progress_reports(rem_progress_writes);

    if (swtch["limited_time"]) timer.allocate_time(available_run_time);

    timer.Start(icyc);

    for (;icyc<MCCYC; ++icyc) {
    	
    	if (debug) std::cout<<"\n\nMAIN LOOP: run cycle "<<icyc<<" of "<<MCCYC<<"\n";
        mc->RunCycle();

        ffh.interaction_potential()->reset_total();

        H.sample(icyc,curTindex);

        if ((icyc+1)%NRT==0) {
            H.writeRTSnapshot();
            PH.population()->SaveSnapshot(snapshot_format, pdbfil,icyc,
                                          curTindex,
                                          ffh.interaction_potential()->value());

            if (ffh.interaction_potential()->value()<emin) {
                emin=ffh.interaction_potential()->value();
                PH.population()->SaveSnapshot(snapshot_format,minpdb,icyc,
                                              curTindex,emin);
            }
        }

        if ((icyc+1)%IAVG==0) {
            H.writeAverages();
            H.writeHistograms();
        }

        if ((icyc+1)==NTHERM && swtch[string("thermalisation")]) {
            prf::cout<<"Requested thermalisation interval of "<<NTHERM
                    <<" cycles has been reached. Starting to collect data "
                    <<"for averages and histograms\n\n";
            H.enableStatistics();
        }

        if ((icyc+1)%ICONF==0) {
            writeConf();
            write_MC_averages();
            timer.record(icyc+1);
            timer.forecast();

            if (timer.recommendation()==Timer::suspend or (icyc+1)==MCCYC) {
                ranh.generator()->saveState(rnboosterfile);
                close_traj_file();
            }
            if (timer.recommendation()==Timer::suspend) {
                prf::cout<<"Suspending job after "<<icyc+1<<" cycles "
                <<"following recommendation from timer.\n";
                break;
            }
        }

        if ((icyc+1)%ncyc_per_T_updt==0) update_T(); // Empty function for Basic MC
    }

    H.writeAverages();
    H.writeHistograms();
    write_MC_averages();
    timer.Stop(icyc+1);
    timer.flush_timing_data(speedfile.c_str());
    prf::cout<<"\n-- Closing STDOUT at the end or suspension of a run --\n\n";
    prf::cerr<<"\n-- Closing STDERR at the end or suspension of a run --\n\n";
    prf::clog<<"\n-- Closing logfile at the end or suspension of a run --\n\n";
    prf::cout.close();
    prf::cerr.close();
    prf::clog.close();
    return 1;
}

void BasicMCRun::auto_track_obs()
{
    ForceField *ff=ffh.interaction_potential();
    H.track(&esum,"rt avg");
    

    for (size_t i=0;i<ff->n_terms();++i) {
    	H.track(ff->term(i),"rt avg");
    }
    if (ffh.interaction_potential()->hasNatCon())
        	H.track(&ego,"rt avg");
    //H.make_obs(InstructionString("RCBin HelixContent default_limits helix"));

    //H.make_obs(InstructionString("RCBin BetaStrandContent default_limits strand"));
}

void BasicMCRun::print_MC_setup()
{
    highlight(progheader.c_str());
    mc->print_setup();
}

void BasicMCRun::printMCconfig()
{
    highlight("System Configuration");
    PH.population()->WriteShort();
    prf::cout<<"Box size "<<AtomCoordinates::boxL()<<" Angstroms\n";
    print_MC_setup();
    mc->updates_hander()->print_updateprobs(prf::cout);
    highlight("Random Number Generator");
    prf::cout<<ranh.generator()->Name()<<" with seed "
            <<ranh.generator()->Seed()<<"\n";
}

#ifndef MAINFUNC

#ifndef PARALLEL
int main(int argc, char *argv[])
{

    BasicMCRun interface;

    if (interface.init(argc,argv)) interface.Run();

    return 0;
}

#else
int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    BasicMCRun interface;
    interface.disable("rank");
    interface.disable("nruns");
    interface.set_rank_nruns(MPI::COMM_WORLD.Get_rank(),
                             MPI::COMM_WORLD.Get_size());
    std::cout<<"starting MPI run:"<<MPI::COMM_WORLD.Get_size()<<"\n";
    if (interface.init(argc,argv)) interface.Run();

    MPI::Finalize();

    return 0;
}

#endif
#endif

int BasicMCRun::init(int argc, char *argv[])
{
    // First check command line
    optn.analyze(argc,argv);
    // Do things that can not wait ...
    if (optn.option_given("rank")) myrank=atoi(optn.option("rank").c_str());
    if (optn.option_given("nruns")) nruns=atoi(optn.option("nruns").c_str());
    if (optn.option_given("debug")) {
    	
    	std::cout<<"options: debug?"<<debug<<"\n";
    }
    optn.set_rank(myrank);
    if (optn.switch_given("help")) {
        for (int i=0;i<optn.n_spare_args();++i) {
            help(optn.spare_args(i));
        }
        if (optn.n_spare_args()==0) show_basic_help();
        return 0;
    }
    // Now we can create the run directory
    init_dir();
    // If an alternative settings file is given ...
    if (optn.option_given("settings_file"))
        optn.settings_file_name(optn.option("settings_file"));
    // Now it is time to get the primary and secondary settings
    optn.get_settings();
    // Add the command line options to the settings file instructions
    optn.add_cmd_line_instructions();
    // Import one combined list of instructions
    std::list<InstructionString> cmds=optn.get_options();
    // Execute instructions with direct meanings inside this class
    parseCommands(cmds);
    // Initialize output streams
    init_streams();
    // Initialize filenames for various output files
    init_filenames();
    icyc0=icyc=0;
    // Sniff and find out if this run starts from a saved state
    int recover_code=recover(RMCCYC);
    if (recover_code==1) {
        swtch["continue"]=true;
        swtch["preliminary_relaxation"]=false;
        prf_time t;
        prf::cout<<"Continue at cycle "<<icyc<<" at "<<t.to_UTC()<<"\n";
        prf::clog<<"Continue at cycle "<<icyc<<" at "<<t.to_UTC()<<"\n";
    } else if (recover_code==0) {
        // Auto-seed and check if there are instructions for the
        //random number generator
        ranh.auto_seed(myrank);
        ranh.parseCommands(cmds,argc,argv);
        // Set up population
        PH.population()->RandomNumberGenerator(ranh.generator());
        PH.parseCommands(cmds,argc,argv);
        if (PH.init_pop()==0) {
            prf::cerr<<"Population initialilsation failed. Population can be set up "
                    <<"using command line arguments or the settings file. \n";
            prf::cerr<<"\nRun "<<argv[0]<<" -h \n\nfor usage information.\n";
            return 0;
        }
        if (PH.init_coords()!=0) swtch["preliminary_relaxation"]=false;
    } else {
        prf::cerr<<"Unsuccessful resume or continue.\n";
        return 0;
    }
    PH.reconstruct();
    // Set up interaction potential
    ffh.parseCommands(cmds,argc,argv);
    
    ffh.init_ff();
    ffh.set_population(PH.population());
    ffh.interaction_potential()->setDebug(debug);
    ffh.interaction_potential()->init();
    ffh.interaction_potential()->evaluate();
    esum.connect(ffh.interaction_potential());
    esum.refresh();
    if (ffh.interaction_potential()->hasNatCon()){
    	std::cout<<"Connect Ego\n";
    	ego.connect(ffh.interaction_potential());
    	ego.refresh();
    }
    
    if (recover_code==1) {
        if (fabs(esum.Value()-en)>1e-8) {
            prf::cerr<<"At cycle "<<icyc<<", \n";
            prf::cerr<<"Energy retrieved from conf file = "<<en
                    <<" Energy calculated anew = "<<esum.Value()<<"\n";
            ffh.interaction_potential()->print_contributions(prf::cerr);
        }
    }
    
    
    // Set up Monte Carlo component
    mc->Connect(PH.population());
    mc->forcefield_handler(&ffh);
    mc->RandomNumberGenerator(ranh.generator());
    mc->updates_hander()->set_n_temps(1);
    mc->parseCommands(cmds,argc,argv);
    int nerr=init_MC();
    Tcur=mc->temperature();
    std::cout << "Tcur="<<Tcur<<"\n";
    if (nerr!=0) {
        prf::cerr<<argv[0]<<"> "<<nerr<<" errors in set up of the "
                <<"Monte Carlo. Can not proceed. \n";
        return 0;
    }
    if (recover_code==1) update_T();
    // Set up measurements
    Logger()(10)<<"Initializing observables and ObsHandler..\n";
    H.setPopulation(PH.population());
    std::cout<<"BasicMCRun.init() ... NTMP="<<NTMP<<"\n";
    H.set_block_props(NTMP,"temperature");
    auto_track_obs();
    H.parseCommands(cmds,argc,argv);
    H.setPrefix(mydir+"/");

    // Some energy terms can depend on things calculated during some measurement
    for (size_t deptrm=0;deptrm<ffh.n_obs_dependent_terms();++deptrm) {
        H.fix_dependencies(ffh.obs_dependent_term(deptrm));
    }

    H.initialize();
    H.writeRTKey();
    writeTemperatures();
    prf::greet();

    return 1;
}

void BasicMCRun::show_basic_help()
{
    prf::cout<<"BasicMCRun: BasicMCRun is a program to do constant temperature\n"
            <<"canonnical Monte Carlo simulations with ProFASi. There are many\n"
            <<"command line options for this program. They are divided into\n"
            <<"different categories. Run the program with the option\n"
            <<"--help category_name(s) to get help on that category. The long\n"
            <<"forms of all command line options can be used as commands in\n"
            <<"a settings file. \n\n"
            <<"Help is available for the following categories ...\n\n"
            <<"Controls: run length, available time, output redirection etc. \n"
            <<"Population: sequence, number of chain, box length etc. \n"
            <<"Energy: Choice of force field \n"
            <<"Random: Choice of random number generator \n"
            <<"MC: Monte Carlo options, temperature, conformational updates etc.\n";
    prf::cout<<"For example, you could write \n\n"
            <<"BasicMCRun --help Population\n\n"
            <<"to see command line options about setting up the population. \n\n";
}

void BasicMCRun::help(std::string qury)
{
    if (qury=="Controls" or qury=="controls") {
        prf::cout<<"\nOptions in the \"Controls\" group ...\n";
        show_help();
    } else if (qury=="Population" or qury=="population") {
        prf::cout<<"\nOptions in the \"Population\" group ...\n";
        PH.show_help();
    } else if (qury=="Energy" or qury=="energy") {
        prf::cout<<"\nOptions in the \"Energy\" group ...\n";
        ffh.show_help();
    } else if (qury=="Random" or qury=="random") {
        prf::cout<<"\nOptions in the \"Random\" group ...\n";
        ranh.show_help();
    } else if (qury=="MC" or qury=="mc") {
        prf::cout<<"\nOptions in the \"MC\" group ...\n";
        mc->show_help();
    }
}


int BasicMCRun::init_dir()
{
    char dirnam[6];
    sprintf(dirnam,"n%d",myrank);
    mydir=string(dirnam);
    CreateDir(mydir);   //check if the directory exists, else create it
    prf::clog.open((mydir+string("/logfile")).c_str(),"a");
    return 1;
}

FILE * BasicMCRun::open_segment(std::string prfx)
{
    prf_time tnow;
    datafl="conf.data_"+tnow.stamp();
    std::string trajfl=prfx+"traj";
    Output trj;
    if (prf_utils::TestFile_r(trajfl)) trj.open(trajfl.c_str(),"a");
    else {
        trj.open(trajfl.c_str(),"w");
        trj<<"#_PROFASI_TRAJECTORY_DESCRIPTION_FILE\n"
                <<"# LAYOUT: SEGMENT_FILE START_CYCLE INTERVAL N_BLOCKS\n";
    }
    trj<<datafl<<"\n";
    trj.close();
    datafl=prfx+datafl;
    FILE *fp=fopen(datafl.c_str(),"w");

    int len=0;
    len+=sizeof(unsigned long)/sizeof(char);
    len+=sizeof(int)/sizeof(char);
    len+=sizeof(double)/sizeof(char);
    int ranbegin=len;
    len+=ranh.generator()->ConfSize();
    int ranend=len;

    for (int j=0;j<PH.population()->NumberOfChains();++j) {
        len+=PH.population()->Chain(j)->ConfSize();
    }

    fprintf(fp,"PROFASI_CONF_HEADER\n");
    fprintf(fp,"version_stamp %s\n",profasi_git_version().c_str());
    fprintf(fp,"conf_length %d\n",len);
    fprintf(fp,"byte_order %s\n", running_on_little_endian() ?
            "little_endian" : "big_endian");

    fprintf(fp, "unsigned long 1 %d %d\n",sizeof(unsigned long)/sizeof(char),
            sizeof(unsigned long)/sizeof(char));
    fprintf(fp, "int 1 %d %d\n",sizeof(int)/sizeof(char),
            sizeof(int)/sizeof(char));
    fprintf(fp, "double 1 %d %d\n",sizeof(double)/sizeof(char),
            sizeof(double)/sizeof(char));
    fprintf(fp,"random_number_start %d\n",ranbegin);
    fprintf(fp,ranh.generator()->ConfSignature().c_str());
    fprintf(fp,"random_number_end %d\n",ranend);
    len=ranend;

    for (int j=0;j<PH.population()->NumberOfChains();++j) {
        fprintf(fp,"peptide_%d_start %d\n",j,len);
        fprintf(fp, PH.population()->Chain(j)->ConfSignature().c_str());
        len+=PH.population()->Chain(j)->ConfSize();
        fprintf(fp,"peptide_%d_end %d\n",j,len);
    }

    fprintf(fp, "box_length %f\n", AtomCoordinates::boxL());
    fprintf(fp, mc->ConfSignature().c_str());

    fprintf(fp,"sizeof_char %d\n",(int) sizeof(char));
    fprintf(fp,"sizeof_int %d\n",(int) sizeof(int));
    fprintf(fp,"sizeof_long %d\n",(int)sizeof(long));
    fprintf(fp,"sizeof_unsigned_long %d\n",(int)sizeof(unsigned long));
    fprintf(fp,"sizeof_float %d\n",(int)sizeof(float));
    fprintf(fp,"sizeof_double %d\n",(int)sizeof(double));
    fprintf(fp,"sizeof_long_double %d\n",(int)sizeof(long double));
    fprintf(fp,"END_PROFASI_CONF_HEADER\n");
    return fp;
}

int BasicMCRun::close_traj_file()
{
    std::string trjfl=mydir+"/traj";
    std::deque<std::string> lines;
    if (prf_utils::get_lines(trjfl,lines)) {
        Output trj(trjfl.c_str());
        for (size_t i=0;i<(lines.size()-1);++i) {
            trj<<lines[i]<<"\n";
        }
        std::deque<std::string> parts;
        prf_utils::split(lines.back(),parts);
        if (!parts.empty()) {
            if ((mydir+"/"+parts[0])==datafl) {
                trj<<parts[0]<<" "<<icyc0<<" "<<ICONF<<" "<<nblk<<"\n";
            }
        }
        trj.close();
    }
    return 1;
}

void BasicMCRun::writeConf()
{
    FILE *fp;
    if (nblk==0) {
        fp=open_segment(mydir+"/");
        icyc0=icyc;
    } else {
        fp=fopen(datafl.c_str(),"a");
    }
    en=esum.Value();
    fwrite(&icyc,sizeof(unsigned long),1,fp);
    fwrite(&curTindex,sizeof(int),1,fp);
    fwrite(&en,sizeof(double),1,fp);
    //if (ffh.interaction_potential()->hasNatCon()){
    //	fwrite(&ego.Value(),sizeof(double),1,fp);
    //}
    ranh.generator()->WriteTo(fp);

    for (int j=0;j<PH.population()->NumberOfChains();++j)
        PH.population()->Chain(j)->WriteConf(fp);

    fclose(fp);
    ++nblk;
}

int BasicMCRun::recover(unsigned rcyc)
{
    prf_traj::Trajectory traj;
    traj.parse(trajfile);
    if (traj.init()==0) return 0;

    unsigned long cycmin=0,cycmax=0;
    get_span(traj,cycmin,cycmax);
    if ((rcyc>cycmax or rcyc<cycmin) and rcyc!=((unsigned)-1)) {
        prf::cerr<<"Resume> The restart point "<<rcyc
                <<" is not within the available range "
                <<cycmin<<" to "<<cycmax<<"\n";
    }
    if (rcyc>cycmax) rcyc=cycmax;
    if (rcyc<cycmin) rcyc=cycmin;
    prf::cerr<<"Resume> Trying to resume at cycle "<<rcyc<<".\n";
    prf_traj::Trajectory::iterator cb=traj.find(rcyc);

    if (cb==traj.end()) {
        prf::cerr<<"Retrieving configuration at cycle "<<rcyc
                <<" from trajectory file "<<trajfile<<" failed.\n";
        return 0;
    } else {
        curTindex=cb->T_index();
        icyc0=cb->MC_time();
        icyc=icyc0+1;
        prf::cout<<"Initializing cycle count to "<<icyc<<"\n";
        en=cb->energy();
    }
    ranh.generator()->set_seed_and_ncalls(cb->random_number_seed(),
                             cb->random_number_ncalls());
    ranh.generator()->retrace(rnboosterfile);
    prf_xml::XML_Node *confmap=traj.xml_map();

    if (confmap!=NULL) {
        if (confmap->child("box_length")!=NULL) {
            PH.boxSize(strtod(confmap->child("box_length")->value().c_str(),
                              NULL),true);
        }
        if (confmap->child("population")!=NULL) {
            PH.population()->assign_sequences(confmap->child("population"));
        }
    }

    if (PH.init_pop()==0) {
        prf::cerr<<"This seems to be a continuation run. Population set up"
                <<"was attempted using stored configuration, and failed. \n"
                <<"This run will exit to draw attention to the problem. To "
                <<"start a fresh run, simply remove or rename the traj files \n"
                <<"n0/traj, n1/traj etc., and try again.\n";
        return -1;
    }
    PH.population()->set_dof(cb->coordinates());

    std::string trajfl=mydir+"/traj";
    if (rcyc!=cycmax and prf_utils::TestFile_r(trajfl)) {
        prf_time tnow;
        std::string bkpid=tnow.stamp();
        rename(trajfl.c_str(),(trajfl+".bkp_"+bkpid).c_str());
        traj.save_up_to(rcyc,trajfl);
        if (prf_utils::TestFile_r(rnboosterfile)) {
            std::string bkpconf=rnboosterfile+".bkp_"+bkpid;
            rename(rnboosterfile.c_str(),bkpconf.c_str());
        }
    }
    return 1;
}

int BasicMCRun::get_span(prf_traj::Trajectory &traj,
                         unsigned long &i1, unsigned long &i2)
{
    i1=traj.min_cycle();
    i2=traj.max_cycle();
    return 1;
}
