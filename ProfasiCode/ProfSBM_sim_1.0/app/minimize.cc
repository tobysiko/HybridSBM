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

#include "minimize.hh"

minimize::minimize()
{
    outputfile="output.xml";
    optn.option("output_file","o",1,
                "(File to save the final conformation)");
    optn.option("using","u",2,"(DOF type and selection)");
    optn.option("log_level","ll",1,
                "(Verbosity level of log messages)");
    optn.option("max_cycles","mxcyc",1,
                "Maximum number of CG cycles");
    sel="*";
    typ="all";
}

minimize::~minimize() {}

int minimize::run()
{
//    PH.population()->SaveSnapshot(1,"start.pdb",0,0,0);
    double fvl=cg.minimize();

    size_t extn=outputfile.find_last_of('.');
    oformat=outputfile.substr(extn+1);
    if (oformat.empty()) oformat="xml";
    int infmt=0;
    if (oformat=="pdb") infmt=1;
    else if (oformat=="xml") infmt=2;
    else if (oformat=="tcn" or oformat=="tconf" or oformat=="tcnf")
        infmt=3;
    else if (oformat=="bcn") infmt=4;

    Logger(10)<<"Writing output file "<<outputfile<<" in "<<oformat<<" format.\n";
    PH.population()->SaveSnapshot(infmt,outputfile,0,0,fvl);
    return 0;
}

void minimize::parseCommands(std::list<InstructionString> &cmds)
{
    for (std::list<InstructionString>::iterator it=cmds.begin();
    it!=cmds.end();++it) {
        parseCommand(*it);
    }
}

void minimize::parseCommand(InstructionString &s)
{
    if (s.head()=="output_file") {
        outputfile=s.tail().str();
    } else if (s.head()=="log_level") {
        int llvl=atoi(s.tail().str().c_str());
        prf::Logger::verbosity=llvl;
    } else if (s.head()=="using") {
        if (s.n_parts()!=3) {
            prf::cerr<<"Incorrect syntax for option using : "<<s.str()<<"\n";
            return;
        }
        typ=s.part(1);
        sel=s.part(2);
    } else if (s.head()=="max_cycles"){
        int ncyc=atoi(s.tail().str().c_str());
        if (ncyc>0) cg.set_max_cyc(ncyc);
    }
}

void minimize::auto_track_obs()
{
    ForceField *ff=ffh.interaction_potential();
    H.track(&esum,"rt avg");

    for (size_t i=0;i<ff->n_terms();++i) H.track(ff->term(i),"rt avg");

    H.make_obs(InstructionString("RCBin HelixContent default_limits helix"));
    H.make_obs(InstructionString("RCBin BetaStrandContent default_limits strand"));
}

int minimize::init(int argc, char *argv[])
{
    prf::Logger::verbosity=20;
    optn.analyze(argc,argv);
    if (optn.n_spare_args()>0) expression=optn.spare_args(0);
    else {
        prf::cerr<<"Did not find expression to minimize in command line.\n";
        return 0;
    }
    if (optn.switch_given("help")) {
        for (int i=0;i<optn.n_spare_args();++i) {
            help(optn.spare_args(i));
        }
        if (optn.n_spare_args()==0) show_basic_help();
        return 0;
    }
    if (optn.option_given("settings_file"))
        optn.settings_file_name(optn.option("settings_file"));
    optn.get_settings();
    optn.add_cmd_line_instructions();
    std::list<InstructionString> cmds=optn.get_options();
    PH.parseCommands(cmds,argc,argv);
    if (PH.init_pop()==0) {
        prf::cerr<<"Population initialilsation failed. Population can be set up "
                    <<"using command line arguments or the settings file. \n";
        prf::cerr<<"\nRun "<<argv[0]<<" -h \n\nfor usage information.\n";
        return 0;
    }
    PH.init_coords();
    PH.reconstruct();
    m.set_population(PH.population());
    ffh.parseCommands(cmds,argc,argv);
    ffh.init_ff();
    ffh.set_population(PH.population());
    ffh.interaction_potential()->init();
    ffh.interaction_potential()->evaluate();
    esum.connect(ffh.interaction_potential());
    esum.refresh();
    ffh.interaction_potential()->print_contributions(prf::cerr);
    H.setPopulation(PH.population());
    H.set_block_props(1,"temperature");
    auto_track_obs();
    H.parseCommands(cmds,argc,argv);
    for (size_t deptrm=0;deptrm<ffh.n_obs_dependent_terms();++deptrm) {
        H.fix_dependencies(ffh.obs_dependent_term(deptrm));
    }

    H.setSwitch("histograms",false);
    H.initialize();

    // Ready to see what needs to be minimized
    std::deque<std::string> terms;
    prf_utils::split_str(expression,'+',terms);
    for (size_t i=0;i<terms.size();++i) {
        std::deque<std::string> parts;
        prf_utils::split_str(terms[i],'*',parts);
        if (parts.size()==1) {
            if (prf_utils::is_number(parts[0])) {
                prf::cerr<<"Ignoring constant "<<parts[0]<<" in minimization.\n";
            } else {
                Observable *o=H.get_obs(parts[0]);
                if (o==NULL) {
                    prf::cerr<<"Unrecognized constant or observable name \""
                            <<parts[0]<<"\" ignored in minimization.\n";
                } else {
                    m.add_obs(1,o);
                }
            }
        } else if (parts.size()==2) {
            bool hasprefactor=false;
            double prefct=0;
            if (prf_utils::is_number(parts[0])) {
                prefct=strtod(parts[0].c_str(),NULL);
                hasprefactor=true;
            } else if (prf_utils::is_number(parts[1])) {
                prefct=strtod(parts[1].c_str(),NULL);
                hasprefactor=true;
            }
            Observable *o=H.get_obs(parts[0]);
            if (o==NULL) o=H.get_obs(parts[1]);

            if (hasprefactor and o!=NULL) {
                m.add_obs(prefct,o);
            }
        } else {
            prf::cerr<<"Could not interpret token "<<terms[i]<<"\n";
        }
    }

    alldofs=PH.population()->dof_map();
    parseCommands(cmds);

    int ich=0,res0=0,res1=0;
    if (sel!="*") {
        std::deque<std::string> parts;
        prf_utils::split_str(sel,'/',parts);
        if (parts.size()>0) {
            ich=atoi(parts[0].c_str());
            if (ich<0 or ich>=PH.population()->num_chains()) return 0;
            res0=0;
            res1=PH.population()->Chain(ich)->numLigands();
        }
        if (parts.size()>1) {
            res0=atoi(parts[1].c_str());
            res1=res0+1;
        }
        if (parts.size()>2) res1=atoi(parts[2].c_str());
        if (res0<0 or res1<0 or res0>=res1 or
            res0>=PH.population()->Chain(ich)->numLigands() or
            res1>=PH.population()->Chain(ich)->numLigands()) return 0;

    }
    for (size_t i=0;i<alldofs.size();++i) {
        DOF_Info dof=alldofs[i];
        bool insel=false;
        if (sel=="*") insel=true;
        else {
            insel=(dof.chain==ich && dof.group>=res0 && dof.group<res1);
        }
        if (!insel) continue;
        if (typ=="all") m.add_dof(alldofs[i]);
        else if (typ=="bb" && dof.dof_kind==backbone_torsion_angle)
            m.add_dof(dof);
        else if (typ=="sc" && dof.dof_kind==sidechain_torsion_angle)
            m.add_dof(dof);
    }
    m.init();
    cg.set_function(&m);
    return 1;
}

void minimize::help(std::string tpc)
{
    prf::cout<<"Help on "<<tpc<<"\n";
}

void minimize::show_basic_help()
{
    prf::cout<<"Some basic help \n";
}

int main(int argc, char *argv[])
{
    minimize interface;
    if (interface.init(argc, argv)) interface.run();
    return 0;
}
