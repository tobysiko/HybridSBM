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

#include "extract_props.hh"
#include <Aux/InstructionString.hh>


extract_props::extract_props()
{
    irtstart=0;irtend=(unsigned long) -1;
    every=1;
    optn.new_switch("get_rt","rt",false,"(whether to reconstruct rt history)");
    optn.new_switch("get_averages","avg",false,"(generage averages file)");
    optn.new_switch("get_pdbs","pdb",false,"(whether to retrive pdb snapshots)");
    optn.new_switch("recenter_all_frames","recenter",false,
                    "(Center each frame around the center of mass (jittery))");
    optn.new_switch("single_file","sf",false,
                    "(whether all snapshots should be written in one file)");
    optn.new_switch("raw", "r",false,
                    "(Read raw binary files instead of traj files)");
    optn.option("start","a",1,"(MC cycle index where reconstruction starts)");
    optn.option("end","z",1,"(MC cycle index where reconstruction ends)");
    optn.option("every","ev",1,"(Use every n'th configuration in file)");
    optn.option("output_prefix","op",1,
                "(prepend for all generated output files)");
    optn.option("number_of_temperatures","nt",1,
                "(number of temperature indexes in the file)");
    optn.option("log_level","ll",1,"(Verbosity of log messages)");

    mytraj="traj";

    NTMP_known=false;
    NTMP=1;
    myprefix="./";
    log_level=3;
}

extract_props::~extract_props() {}

void extract_props::show_basic_help()
{
    prf::cout<<"\nThis program interprets trajectory information from "
            <<"PROFASI's traj files and regenerates requested information, "
            <<"such as PDB snapshots in a given interval, the rt file in a"
            <<"certain range, averages etc.\n\nAvailable options:\n\n";
    optn.write_available();
    prf::cout<<"Examples:\n\n"
            <<"extract_props -pdb -a 39999 -z 59999 n0/traj -op tmp_ \n\n"
            <<"The above takes all configurations between MC cycles 39999 and "
            <<"59999 and saves a PDB file for each. The DOF information at "
            <<"those cycles, is stored in some 'conf.data...' file referred "
            <<"to by the traj file. The PDB files are named with a prefix "
            <<"tmp_ : tmp_frame_1.pdb, tmp_frame_2.pdb etc. \n\n"
            <<"extract_props -rt -avg -a 39999 -z 59999 n0/traj -op tmp_ \n\n"
            <<"The above takes the same range of cycles as in the previous "
            <<"example, and generates rt and averages files.\n\n";
    prf::cout<<"For more information on more options available for this "
            <<"program, consult the documentation.\n";
}

int extract_props::my_init(int argc, char *argv[])
{
    // First check command line
    optn.analyze(argc,argv);
    if (argc==1 or optn.switch_given("help")) {
        show_basic_help();
        return 0;
    }
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
    for (std::list<InstructionString>::iterator it=cmds.begin();
    it!=cmds.end();++it) {
        parseCommand(*it);
    }
    prf::Logger::verbosity=log_level;

    if (optn.n_spare_args()<1) {
        prf::cerr<<"No trajectory file given.\n";
        return 0;
    }

    if (!optn.switch_given("r")) {
        std::string inputconf=optn.spare_args(0);

        int errcode=traj.parse(inputconf);
        if (errcode==0) {
            prf::cerr<<"No valid segments in trajectory.\n";
            return 0;
        } else if (errcode==-1) {
            prf::cerr<<"Failed to interpret file "<<inputconf<<".\n";
            prf::cerr<<"If you were trying to extract from a file containing "
                    <<"binary configuration data, use the option \"--raw\" with "
                    <<"the appropriate syntax.\n";
            return 0;
        }
    } else {
        std::deque<std::string> inputconf;
        for (int i=0;i<optn.n_spare_args();++i)
            inputconf.push_back(optn.spare_args(i));
        traj.append_list(".",inputconf);
    }
    if (traj.init()==0) return 0;

    prf_xml::XML_Node *confmap=traj.xml_map();

    if (confmap!=NULL) {
        if (confmap->child("box_length")!=NULL) {
            PH.boxSize(strtod(confmap->child("box_length")->value().c_str(),
                              NULL),true);
        }
        if (confmap->child("population")!=NULL) {
            PH.population()->assign_sequences(confmap->child("population"));
        }
        if (confmap->child("num_temperatures")!=NULL) {
            NTMP=atoi(confmap->child("num_temperatures")->value().c_str());
            NTMP_known=true;
        }
    }

    if (not NTMP_known) {
        std::string tinfo=optn.spare_args(0);
        unsigned loc=tinfo.find_last_of('/');
        if (loc>=tinfo.size()) tinfo="conf.info";
        else tinfo=tinfo.substr(0,loc)+"/temperature.info";
        std::deque<std::string> tlines;
        prf_utils::get_lines(tinfo,tlines);
        if (tlines.size()>1) {
            NTMP=tlines.size()-1;
            NTMP_known=true;
        }
    }

    if ((optn.switch_given("rt") or
         optn.switch_given("avg") or
         optn.switch_given("his")) and
        not NTMP_known) {
        prf::cerr<<"For what you are trying to do, the number of temperature "
                <<"indexes to be expected in the trajectory files needs to be "
                <<"known. I can't figure it out. Please rerun with the option "
                <<"-nt or --number_of_temperatures\n";
        return 0;
    }

    if (PH.init_pop()==0) {
        prf::cerr<<"Population initialilsation failed. Probably the binary "
                <<"configuration file and the optional accompanying info file "
                <<"did not have enough information to build the population. ";
        prf::cerr<<"Try generate a conf.info file using the program "
                <<"generate_conf_info, and pass that to this program using the "
                <<"--info_file option\n";
        return 0;
    }
    PH.init_coords();

    // Set up interaction potential
    ffh.parseCommands(cmds,argc,argv);
    ffh.init_ff();
    ffh.set_population(PH.population());
    ffh.interaction_potential()->init();
    esum.connect(ffh.interaction_potential());

    unsigned long cstart=traj.min_cycle(),cend=traj.max_cycle();

    Logger(5)<<"The given trajectory spans from "<<cstart<<" to "<<cend<<"\n";
    if ((!optn.option_given("a")) or irtstart<cstart) {
        irtstart=cstart;
        if (optn.option_given("a"))
            prf::cerr<<"Using the minimum available cycle "<<cstart
                <<" as start cycle.\n";
    }
    if ((!optn.option_given("z")) or irtend>cend) {
        irtend=cend;
        if (optn.option_given("z"))
            prf::cerr<<"Using the maximum available cycle "<<cend
                <<" as end cycle.\n";
    }

    cb=traj.find(irtstart);
    if (cb==traj.end()) {
        prf::cerr<<"Cycle "<<irtstart<<" couldn't be found.\n";
        return 0;
    }

    PH.population()->set_dof(cb->coordinates());
    PH.reconstruct();
    ffh.interaction_potential()->evaluate();
    esum.refresh();

    // Set up measurements
    Logger()(10)<<"Initializing observables and ObsHandler..\n";
    H.setPopulation(PH.population());
    H.unsetSwitch("histograms");
    H.set_block_props(NTMP,"temperature");
    auto_track_obs();
    H.parseCommands(cmds,argc,argv);
    H.initialize();
    H.setPrefix(myprefix);

    if (optn.switch_given("rt")) H.writeRTKey();

    return 1;
}

void extract_props::auto_track_obs()
{
    ForceField *ff=ffh.interaction_potential();
    H.track(&esum,"rt avg");

    for (size_t i=0;i<ff->n_terms();++i) H.track(ff->term(i),"rt avg");

    H.make_obs(InstructionString("RCBin HelixContent default_limits helix"));

    H.make_obs(InstructionString("RCBin BetaStrandContent default_limits strand"));
}

int extract_props::parseCommand(InstructionString s)
{
    if (s.head()=="output_prefix") myprefix=s.tail().str();
    else if (s.head()=="log_level") log_level=atoi(s.tail().str().c_str());
    else if (s.head()=="start") irtstart=strtoul(s.tail().str().c_str(),NULL,10);
    else if (s.head()=="end") irtend=strtoul(s.tail().str().c_str(),NULL,10);
    else if (s.head()=="every") every=atoi(s.tail().str().c_str());
    else if (s.head()=="number_of_temperatures") {
        NTMP=atoi(s.tail().str().c_str());
        NTMP_known=true;
    }
    return 1;
}

int extract_props::BrowseConfs()
{
    Logger(5)<<"Requested starting MC time = "<<irtstart<<"\n";

    if (fabs(esum.Value()-cb->energy())>1e-8) {
        prf::cerr<<"At cycle "<<cb->MC_time()<<", \n";
        prf::cerr<<"Energy retrieved from file = "<<cb->energy()
                <<" Energy calculated anew = "<<esum.Value()<<"\n";
        ffh.interaction_potential()->print_contributions(prf::cerr);
    }

    int ndata=0;
    FILE *animation=NULL;
    unsigned long i=cb->MC_time();
    bool recenter = optn.switch_given("recenter");

    if (optn.switch_given("pdb") && optn.switch_given("sf")) {
        animation=fopen((myprefix+"frames.pdb").c_str(),"w");
        PH.population()->writePDBHeader(animation,i,cb->T_index(),cb->energy());
    }

    while (i<=irtend) {
        if (ndata%500==0)
            prf::cerr<<"i= "<<i<<"\n";
        ++ndata;

        if (optn.switch_given("rt") or optn.switch_given("avg"))
            H.sample(i,cb->T_index());

        if (optn.switch_given("rt")) H.writeRTSnapshot();

        if (optn.switch_given("pdb")) {
            if (recenter) AtomCoordinates::center();
            PH.population()->EnforceBC();
            AtomCoordinates::update(0,PH.population()->NumberOfAtoms());

            if (!optn.switch_given("sf")) {
                char newpdb[120];
                sprintf(newpdb,"frame_%d.pdb",ndata);
                PH.population()->SaveSnapshot(1,(myprefix+newpdb),i,
                                              cb->T_index(),cb->energy());
            } else {
                fprintf(animation,"MODEL %d\n",ndata);
                PH.population()->WritePDB(animation);
                fprintf(animation,"ENDMDL\n");
            }
        }

        cb+=every;
        cb.refresh();
        if (cb.good()) {
            PH.population()->set_dof(cb->coordinates());
            PH.population()->Reconstruct();
            PH.population()->EnforceBC();
            AtomCoordinates::update(0,PH.population()->NumberOfAtoms());
            i=cb->MC_time();
            ffh.interaction_potential()->evaluate();
            esum.refresh();
        } else break;
    }

    if (optn.switch_given("pdb") && optn.switch_given("sf")) {
        fclose(animation);
    }

//    if (optn.option_given("his")) H.writeHistograms();

    if (optn.switch_given("avg")) H.writeAverages();

    return 1;
}

int main(int argc, char *argv[])
{
    extract_props rtm;

    if (rtm.my_init(argc,argv)) rtm.BrowseConfs();

    return 0;
}
