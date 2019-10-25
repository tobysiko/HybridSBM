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
#include <Aux/fileutils.hh>
#include <Observables/ProteinRMSD.hh>

using namespace prf;

using namespace prf_utils;
using std::string;
using std::vector;
using std::list;
using std::deque;

void write_help(string);
int main(int argc, char *argv[])
{
    ProgArgs par;
    par.option("using","u",1, "(filter for atom catetories)");
    par.option("r_select","rs",1,"(selection for the right hand structure(s))");
    par.option("no_sequence_allign","no",0,
               "(user takes all responsibility for sequence matching)");
    par.option("log_level","ll",1,"(verbosity of log messages)");
    par.option("logfile","l",1);
    par.option("verbose","v",0,"(use it to have more text in standard output)");
    par.option("help","h",0);

    par.analyze(argc,argv);

    if (par.option_given("h")) {write_help(argv[0]);return 0;}

    if (argc==1) {
        prf::cout<<"MimiqA: Minimaler mittlerer quadratischer Abstand\n\n"
        <<"This program calculates root mean square deviation of one "
        <<"or more PDB files with respect to a reference file.\n\n";
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] file1.pdb file2.pdb\n\n";
        prf::cout<<"Options could be any number and combination of ...\n";
        par.write_available();
        prf::cout<<argv[0]<<" file1.pdb file2.pdb -u \"+BB+CB\"\n"
                <<"\nIn the above, file1.pdb is used as the reference, and "
                <<"file2.pdb as the comparison file. There can be any number "
                <<"of comparison files, up to the maximum number of allowed "
                <<"command line arguments on your operating system. \n\n"
        <<argv[0]<<" -u \"+HV-@S\" file1.pdb:2:A,2,23 file2.pdb:1:B,32,53 \n"
        <<argv[0]<<" file1.pdb::A file2.pdb::A\n"
        <<argv[0]<<" -u \"+BB\" nmr.pdb:4:A myfile*.pdb\n\n";
        prf::cout<<"Run \""<<argv[0]
        <<" --help\" for more detailed usage information\n";
        return 0;
    }

    // Some initialisation
    string lgfl=".profasi_logfile";

    Logger blog;

    int loglvl=10;

    if (par.option_given("l")) lgfl=par.option("l");

    prf::clog.open(lgfl.c_str());

    if (par.option_given("ll")) loglvl=atoi(par.option("ll").c_str());

    prf::Logger::verbosity=loglvl;

    prf::Groups::initGroups();

    string rmsdfilter;

    deque<string> files;

    rmsdfilter="+HV";

    if (par.option_given("u")) rmsdfilter=par.option("u");

    if (par.n_spare_args()<2) {
        prf::cerr<<"At least two input pdb files must be given.\n";
        return 1;
    }

    for (int i=0;i<par.n_spare_args();++i) {
        size_t icol=par.spare_args(i).find(':');
        string flnm=par.spare_args(i);

        if ((unsigned) i<par.spare_args(i).size()) flnm=flnm.substr(0,icol);

        if (TestFile_r(flnm.c_str())) files.push_back(par.spare_args(i));
    }

    ProteinRMSD handler;

    if (par.option_given("rs")) {
        handler.set_default_sel2(par.option("rs"));
    }

    handler.set_struc1(files[0]);

    handler.live_first_struc(); //Necessary for possible multiple comparisons

    handler.filter(rmsdfilter);

    if (par.option_given("no")) handler.no_seq_align();

    for (size_t i=1;i<files.size();++i) {
        handler.delete_matrix2();
        handler.set_struc2(files[i]);

        if (handler.init()==0) continue;

        handler.refresh();

        if (par.option_given("v")) prf::cout<<"RMSD("<<files[0]<<" versus "
            <<files[i]<<") = ";

        prf::cout<<handler.Value()<<"\n";
    }

    return 0;
}

void write_help(string prognm)
{
    prf::cout<<"MimiqA: Minimaler mittlerer quadratischer Abstand \n\n";
    prf::cout<<"This program calculates the Root Mean Square Deviation\n"
    <<"between two pdb files. The basic usage is \n\n"
    <<prognm<<" file1 file2\n\n"
    <<"This suffices if the files contain only one model and one\n"
    <<"protein chain each, and the entire range of amino acids\n"
    <<"is to be used. This also assumes that you are interested\n"
    <<"in the RMSD with all heavy (non-hydrogen) atoms.\n"
    <<"If you wish to perform RMSD over backbone atoms you would \n"
    <<"pass a \"using\" option...\n\n"
    <<prognm<<" file1 file2 -u \"+BB\"\n\n"
    <<"When an explicit \"using\" option is given, the default \n"
    <<"filter which chooses heavy atoms is removed and only the \n"
    <<"filter passed explicitly is used. To get RMSD with backbone\n"
    <<"and Cbeta atoms, you would write,\n\n"
    <<prognm<<" file1 file2 -u \"+BB+CB\"\n\n"
    <<"All heavy atoms except those on the backbone ?\n\n"
    <<prognm<<" file1 file2 -u \"+HV-BB\"\n\n"
    <<"The full set of filter and other selection rules is given in "
    <<"the MimiqA page in PROFASI documentation.\n\n"
    <<"To compare one given pdb file with a whole set of files, use\n\n"
    <<prognm<<" nmr.pdb file1 file2 file3 file4 ...\n\n"
    <<"where the dots above are only for effect! \n\n"
    <<"A PDB file may have several chains. To use chain A in one file\n"
    <<"and chain B in the other, you would write,\n\n"
    <<prognm<<" file1::A file2::B -u \"+BB\"\n\n"
    <<"In case there are many chains, and none is selected, all will\n"
    <<"be used. By default, the full range of amino acids in every \n"
    <<"chain are used. If one chain is bigger than the other, the \n"
    <<"smaller one will be slid over the bigger one to find the \n"
    <<"position where the sequence of the smaller one matches.\n"
    <<"This sequence allignment is always done, if the chains differ\n"
    <<"in size or in sequence. For chains with different sequences\n"
    <<"only the parts with overlapping sequence in the best possible\n"
    <<"allignment are used.\n\n"
    <<"For a finer control over the range of residues over which to\n"
    <<"perform calculations, one provides the range along with the \n"
    <<"file name...\n\n"
    <<prognm<<" file1::A,2,33 file2::B,32,63 -u \"+BB\"\n\n"
    <<"It is possible to have discontinuous ranges...\n\n"
    <<prognm<<" file1::A,2,13:A,56,61 file2::A,2,13:A,56,61 -u \"+BB\"\n\n"
    <<"In the above case, residues 2 through 13 and 56 through 61\n"
    <<"in the files will be used. When using disjoint ranges, the \n"
    <<"sequence allignment provided here should not be trusted.\n"
    <<"Instead, the user should carefully choose the segemnts\n"
    <<"from both files.\n\n"
    <<"When a file has several models, one particular model can be \n"
    <<"optionally selected (default: first model)\n\n"
    <<prognm<<" file1:3:A file2 \n\n"
    <<"To print more text on the standard output, pass the -v option\n"
    <<"To increase the verbosity of log messages written in the logfile\n"
    <<"use the --log-level option. Higher values passed to --log-level\n"
    <<"increase verbosity. The log messages are written in the file \n"
    <<"\".profasi_logfile\". This can be changed using the --logfile\n"
    <<"option.\n";
}

/**
\page mimiqa MimiqA : A swiss army knife for RMSD calculations
MimiqA is a command line RMSD calculator. The name stands for "Minimaler mittlerer
 quadratischer Abstand", which is German for "minimal mean square distance". The
 name does not contain the "square root", but the program returns values with the
 square root. MimiqA is implemented using PROFASI version 1.1 PDB Handling
 module. It is not crucial for the other parts of PROFASI to function, and is not
 as such a "Simulation" program. But it is an application with a lot of potential
 use in analyzing protein structures, that extend beyond the intended use of
 PROFASI simulation routines.

Usage: \n \n
\$ mimiqa [OPTIONS] pdbfile_1:model:selections_in_file_1 pdbfile_2:model:selections_in_file_2
\n \n

The program only supports a few options. \n

\li \b --using or \b -u Takes 1 argument, a filter for atom catetories.
\li \b --log_level or \b -ll, Takes 1 argument to set the verbosity of log messages
\li \b --logfile or \b -l, Takes 1 argument, filename, where log messages are stored
\li \b --no_seq_alignment or \b -no Use to prevent attempts to align sequences.
\li \b --verbose or \b -v Use it to have more text in standard output
\li \b --help or \b -h Use it to get a detailed help message on the command line.

While on surface there are only a few options, the program is very flexible. The selections method hinted above, along with the filters passed through the "-u" option allow you to calculate almost any kind of RMSD you might need. We illustrate with a series of examples... \n \n

\section examples
\$ mimiqa file1 file2 \n \n

This suffices if the files contain only one model and one protein chain each, and the entire range of amino acids is to be used. This also assumes that you are interested in the RMSD with all heavy (non-hydrogen) atoms.

If you wish to perform RMSD over backbone atoms you would pass a "using" option...\n \n
\$ mimiqa file1 file2 -u "+BB" \n \n

When an explicit "using" option is given, the default filter which chooses heavy atoms is removed and only the filter passed explicitly is used.

To get RMSD with backbone and Cbeta atoms, you would write, \n \n
\$ mimiqa file1 file2 -u "+BB+CB" \n \n
All heavy atoms except those on the backbone ?\n \n
\$ mimiqa file1 file2 -u "+HV-BB" \n \n

All carbon atoms ?\n \n
\$ mimiqa file1 file2 -u "+@C" \n \n

All backbone atoms but excluding the parts in Proline residues ?\n \n
\$ mimiqa file1 file2 -u "+BB-%PRO" \n \n

You get the idea. The full set of filter and other selection rules is given in \ref prf_sel_fils . \n \n

To compare one given pdb file with a whole set of files, use \n \n
\$ mimiqa nmr.pdb file1 file2 file3 file4 ... \n \n
where the dots above are only for effect! \n \n

A PDB file may have several chains. To use chain A in one file
and chain B in the other, you would write,\n \n
\$ mimiqa file1::A file2::B -u "+BB" \n \n

In case there are many chains, and none is selected, \e ALL will be used.

By default, the full range of amino acids in every chain are used. If one chain is bigger than the other, the smaller one will be slid over the bigger one to find the position where the sequence of the smaller one matches. This sequence allignment is done, when the chains differ in size or in sequence. For chains with different sequences only the parts with overlapping sequence in the best possible allignment are used.\n \n

For a finer control over the range of residues over which to perform calculations, one provides the range along with the file name...\n \n
\$ mimiqa file1::A,2,33 file2::B,32,63 -u "+BB" \n \n

It is possible to have discontinuous ranges...\n \n
\$ mimiqa file1::A,2,13:A,56,61 file2::A,2,13:A,56,61 -u "+BB" \n \n
In the above case, residues 2 through 13 and 56 through 61 in the files will be used. When using disjoint ranges, the sequence allignment provided here should not be trusted. Instead, the user should carefully choose the segemnts from both files.\n \n

When a file has several models, one particular model can be optionally selected (default: first model)\n \n
\$ mimiqa file1:3:A file2 \n \n

Sometimes, it is necessary to compare structurally similar areas in two chains, which have different sequences. Normally, mimiqa would complain about the sequence mismatch, and attempt a sequence alignment. All mismatched parts are ignored. But if we want to use the mismatched parts, say, to calculate backbone RMSD, consisting of atoms common to all residue types, it is convenient to be able to switch-off sequence alignment. The user then takes care of providing an exact range for each structure, and asks PROFASI to not bother with alignment. This is done by passing an option "-no". \n \n
\$ mimiqa -u "+BB" -no mol1.pdb::A,2,25 mol1_mutant.pdb::A,2,25 \n \n
This finds the backbone RMSD in the specified region, without looking at the sequences.

\section notes Additional notes
<ul>
<li> The calculation proceeds in the following order:
<ol>
<li> The PDB files are parsed.</li>
<li> A certain range of residues is selected and the atom records for all atoms in the selected residues are exported </li>
<li> A filter rule is applied on the selected atoms, to reject atoms except those conforming to some criterion.</li>
 <li> The selected atoms of the two structures are compared and the proper correspondence for each atom determined. If one atom does not find a corresponding atom in the other file, it is eliminated.</li>
 <li> The minimal root mean square deviation, the root mean square deviation minimized over all possible relative rotations and translations of the structures,  is calculated between the coordinates of the corresponding atoms using closed form algebraic expressions.</li>
</ol>
</li>

<li> It does not matter if one of the structure files has a few residues somewhere with a few atoms missing. The corresponding atoms in the other structure will be omitted from the comparison.</li>

<li> It does not matter if the PDB files have the atoms ordered differently. The correspondence of each atom is determined during the run. </li>

<li> The log file contains an atom by atom comparison list. So, you can check what is being compared for the RMSD. </li>
</ul>
*/
