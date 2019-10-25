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

#include "extract_snapshot.hh"
#include <Aux/InstructionString.hh>
#include <Aux/profasi_io.hh>
#include <Aux/prf_time.hh>
#include <Aux/profasi_version.hh>

extract_snapshot::extract_snapshot()
{
    optn.settings_use_type(0);
    optn.option("output_file","o",1);
    optn.option("cycle_number", "c",1);
    optn.option("output_format","f",1,"{pdb/xml/binary/textconf}");
    optn.new_switch("stdin", "i", false, "(Take cycle numbers from stdin)");
    optn.option("box_length","box",1,"(Explicitly adjust enclosing box)");
    optn.new_switch("raw","r",false,
                    "(Use binary conf files directly (Inefficient!))");
    std::string formatlist[]={"xml","pdb","bcn","tcn","tconf"};
    formats.clear();
    for (int i=0;i<5;++i) formats.push_back(formatlist[i]);

    outfile="output.xml";
    oformat="xml";
    cycno=0;
    prf::Logger::verbosity=2;
}

extract_snapshot::~extract_snapshot() {}

int extract_snapshot::my_init(int argc, char *argv[])
{
    optn.analyze(argc,argv);
    if (optn.option_given("h") or argc==1) {
        show_help();
        return 0;
    }

    if (!act_on_CL_options()) return 0;

    if (traj.init()) {
        prf_xml::XML_Node *confmap=traj.xml_map();
        double bx=100;
        if (confmap->child("box_length")!=NULL) {
            bx=strtod(confmap->child("box_length")->value().c_str(),NULL);
        } else if (optn.option_given("box")) {
            bx=strtod(optn.option("box").c_str(),NULL);
        }

        PH.boxSize(bx,true);
        PH.population()->assign_sequences(confmap->child("population"));
        if (PH.init_pop()) return 1;
    }

    prf::cerr<<"Initialisation failed\n";

    return 0;
}

void extract_snapshot::show_help()
{
    prf::cout<<"This is a program to extract the instantaneous state or "
            <<"conformation of a system of protein chains, from the trajectory"
            <<" information written by PROFASI simulation programs. Data is "
            <<"stored in a series of binary \"conf...\" files, and metadata "
            <<"about various segments of a run is stored in the trajectory "
            <<"file \"traj\". You should use the \"traj\" files "
            <<"for structure extraction.\n\n";
    prf::cout <<"Example usage: \n\n"
            <<"extract_snapshot -c 4398999 n21/traj -o structure.pdb \n\n"
            <<"extract_snapshot -c 4398999 n21/traj -o structure.xml \n\n"
            <<"extract_snapshot --raw -c 4398999 n21/conf.bkp0 n21/conf.bkp1 "
            <<"n21/conf.bkp2 -o structure.xml\n\n";
    prf::cout<<"The option \"--cycle_number\" or \"-c\" is used to specify "
            <<"the MC cycle number at which information should be retrieved.\n";
    prf::cout<<"The format of the output file is determined by the file "
            <<"extension : .pdb=>pdb, .xml=>xml, .bcn=>binary, "
            <<".tcnf or .tconf => textconf\n";
    prf::cout<<"\n\nThe switch \"raw\" can be used to work in absence of a "
            <<"metadata file, i.e., using directly the binary files, somewhat"
            <<" like in older versions of PROFASI. The difference is that "
            <<"multiple binary conf files can be used at once as in the example"
            <<" above. The raw switch appears convenient, but it is inefficient,"
            <<" and should be avoide when possible. If your data was generated "
            <<"using a version of PROFASI which did not create the metadata "
            <<"trajectory files during the run, you can create them using the "
            <<"program \"create_traj_file\".\n\n";
    prf::cout <<"Here is a summary of supported options ...\n";
    optn.write_available();
}

int extract_snapshot::act_on_CL_options()
{
    if (optn.option_given("o")) outfile=optn.option("o");
    size_t extn=outfile.find_last_of('.');
    oformat=outfile.substr(extn+1);

    if (optn.option_given("f")) oformat=optn.option("f");
    if (oformat=="binary") oformat="bcn";
    if (oformat=="textconf" or oformat=="text_conf") oformat="tcn";

    if (std::find(formats.begin(),formats.end(),oformat)==formats.end())
        oformat="xml";

    if (oformat==outfile.substr(extn+1)) outfile=outfile.substr(0,extn);

    if (optn.option_given("c"))
        cycno=strtoul(optn.option("c").c_str(),NULL,10);

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
    return 1;
}


int extract_snapshot::execute()
{
    std::vector<unsigned long> listcyc;

    if (optn.option_given("c")) listcyc.push_back(cycno);

    if (optn.switch_given("i")) {
        while (std::cin>>cycno) {
            listcyc.push_back(cycno);
        }

        std::sort(listcyc.begin(),listcyc.end());
    }

    std::string ofile=outfile+"."+oformat;

    for (size_t i=0;i<listcyc.size();++i) {
        cycno=listcyc[i];
        prf_traj::Trajectory::iterator cb=traj.find(cycno);
        if (cb!=traj.end()) {
            if (optn.switch_given("i")) {
                char tag[20];
                sprintf(tag,".%lu",listcyc[i]);
                ofile=outfile+std::string(tag)+"."+oformat;
            }
            PH.population()->set_dof(cb->coordinates());
            PH.reconstruct();

            Logger(10)<<"Writing file "<<ofile<<" in "<<oformat<<" format.\n";

            int infmt=0;
            if (oformat=="pdb") infmt=1;
            else if (oformat=="xml") infmt=2;
            else if (oformat=="tcn" or oformat=="tconf" or oformat=="tcnf")
                infmt=3;
            else if (oformat=="bcn") infmt=4;

            PH.population()->SaveSnapshot(infmt,ofile,cb->MC_time(),
                           cb->T_index(),cb->energy());
        }
    }

    return 0;
}

int main(int argc, char *argv[])
{

    extract_snapshot interface;

    if (interface.my_init(argc, argv)) return interface.execute();

    return 1;
}
