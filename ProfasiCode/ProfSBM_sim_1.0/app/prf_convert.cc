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

#include <Elements/PopulationHandler.hh>
#include <Aux/fileutils.hh>
#include <Aux/profasi_version.hh>
#include "Aux/prf_time.hh"
#include <fstream>

using namespace prf;

using namespace prf_utils;

int main(int argc, char *argv[])
{
    ProgArgs prog;
    prog.settings_use_type(2);
    PopulationHandler ph;
    prog.option("output_file","o",1);
    prog.option("log_level","ll",1,"(verbosity level of log messages)");

    prog.init_options(argc,argv);

    if (prog.n_spare_args()==0) {
        prf::cout<<"This program converts between the following file formats:\n"
        <<"(1) PROFASI XML structure format\n"
        <<"(2) PDB format\n"
        <<"(3) PROFASI text configuration format\n\n"
        <<"Usage: "<<argv[0]<<" [OPTIONS] input_file -o output_file\n"
        <<"or, "<<argv[0]<<" [OPTIONS] input_file output_file\n"
        <<"Where OPTIONS could be one or more of ...\n";
        prog.write_available();
        prf::cout<<"Options relating to population set up ...\n";
        ph.show_help();
        prf::cout<<"Note that if the input is in the text configuration format"
        <<", it is necessary to provide additional information so that the "
        <<"numbers in the text configuration file can be interpreted correctly."
        <<"PROFASI text configuration format only stores the degrees of "
        <<"freedom as bare numbers. They only make sense when one knows what "
        <<"system they refer to. So, if that is the input format, you have "
        <<"to use the -st option to specify a settings file, or use the population"
        <<" set up related commands of the setting file as command line options. "
        <<"For instance, you could use a --add_chain option on the command line"
        <<"to tell the converter what system to use to interpret the tconf files.\n"
        <<"If you use the --add_chain option, you have to put the arguments to this "
        <<"option in quotes (see the third example below). All such options are "
        <<"ignored if the structure input format is not tconf.\n\n"
        <<"The format for input and output is inferred from the file extensions. \n\n";
        prf::cout<<"EXAMPLES:\n\n"
        <<argv[0]<<" abc.xml abc.pdb\n\n"
        <<argv[0]<<" abc.tconf -st ../settings.cnf -o abc.xml \n\n"
        <<argv[0]<<" abc.tconf --add_chain 1 \"< *GEWTY DDATKT FTVTE*>\" -o abc.pdb\n\n";
        return 0;
    }

    std::string ofile("output.pdb"), ifile(prog.spare_args(0));

    int loglevel=10;

    if (prog.option_given("ll")) loglevel=atoi(prog.option("o").c_str());

    if (prog.option_given("o")) ofile=prog.option("o");
    else if (prog.n_spare_args()==2) ofile=prog.spare_args(1);

    prf::Logger::verbosity=loglevel;

    std::string iformat,oformat,ifilnm,ofilnm,isel,osel;
    prf_utils::analyze_filename(ifile,ifilnm,iformat,isel);
    prf_utils::analyze_filename(ofile,ofilnm,oformat,osel);

    std::list<InstructionString> cmds;

    if (iformat=="tconf") {
        cmds=prog.get_options();
        if (ifile.substr(0,7)!=std::string("file://")) ifile="file://"+ifile;
        cmds.push_back(InstructionString("init_config "+ifile));
    } else if (iformat=="xml" or iformat=="pdb") {
        cmds.push_back(InstructionString("set_population "+ifile));
    } else {
        prf::cerr<<"Unknown or unsupported input format.\n";
        return 1;
    }
    if (not (oformat=="xml" or oformat=="pdb" or oformat=="tconf" or
        oformat=="tcnf")) {
        prf::cerr<<"Unknown output format.\n";
        return 1;
    }
    ph.parseCommands(cmds,argc,argv);
    if (iformat=="tconf") {
        ph.init_pop();
        ph.init_coords();
        ph.reconstruct();
    }
    Population *p=ph.population();

    FILE *fp=fopen(ofile.c_str(),"w");
    if (oformat=="pdb") {
        unsigned int remn=1;
        fprintf(fp,"REMARK%4u PROGRAM : PROFASI VERSION %s\n",remn++,
        profasi_version().c_str());
        p->writeSequenceInfo(fp);
        p->WritePDB(fp);
    } else if (oformat=="xml") {
        fprintf(fp,"<?xml version=\"1.0\"?>\n<structure>\n");
        fprintf(fp,"<profasi_version>%s</profasi_version>\n",
                profasi_version().c_str());
        fprintf(fp,"<creation_time>\nUTC %s</creation_time>\n",
            prf_time().to_UTC().c_str());
        fprintf(fp,"<box_length>%.16f</box_length>\n",AtomCoordinates::boxL());
        fprintf(fp,"<remark>Converted from %s</remark>\n",ifile.c_str());
        p->Write_XML(fp);
        fprintf(fp,"</structure>\n");
    } else if (oformat=="tconf" or oformat=="tcnf") {
        p->WriteConf_text(fp);
    }
    fclose(fp);

    return 0;
}

/**
  * \page prf_convert Converting between structure file formats
  *
  * The utility \tt prf_convert can convert between the ProFASi XML structure
  * format, ProFASi text configuration format and the PDB format.
  *
  * Examples:
  * <tt>prf_convert abc.xml abc.pdb</tt>\n
  * <tt>prf_convert abc.pdb abc.xml</tt>\n
  * <tt>prf_convert abc.tconf -o abc.pdb --add_chain 10 "ACE*KLVFFAE*NH2"</tt>\n
  * <tt>prf_convert abc.tconf -o abc.xml --settings_file my_settings.cnf</tt>\n
  *
  * \note If neither input nor output format is the tconf format, the input and
  * output filenames are just the first and second command line arguments
  * given. If a tconf file is involved, more arguments are needed, and the
  * user has to specify which of these arguments is the output file name. Also,
  * since the tconf format does not store sequence information, instructions
  * about population set up have to be given either as command line arguments
  * or in a settings file as demonstrated above.
  */
