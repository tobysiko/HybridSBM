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
#include <Elements/PopulationHandler.hh>
#include <Energy/FFHandler.hh>
#include <Observables/ObsHandler.hh>

#include <sstream>

using namespace prf;
using namespace prf_utils;
using namespace prf_pdb_vars;
using namespace std;

class torsion_transforms;

int main(int argc, char *argv[])
{
    ProgArgs par;
    par.settings_use_type(2);
    PopulationHandler ph;
    FFHandler fh;
    par.option("log_level","ll",1);
    par.option("transformations_file","tf",1,"(transform angles with a set of rules)");
    par.init_options(argc,argv);

    if (argc==1) {
        prf::cout<<"This program reads a protein or population structure and prints the "
                <<"ProFASi energy terms. The structure can be given as a PDB file, a "
                <<"ProFASi XML or text structure configuration file.\n\n";
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] structure_file\n\n";
        prf::cout<<"Options could be any of the following ...\n";
        par.write_available();
        prf::cout<<"Options relating to the force field...\n";
        fh.show_help();
        prf::cout<<"Options relating to the population set up\n";
        ph.show_help();
        prf::cout<<"Examples: \n\n";
        prf::cout<<argv[0]<<" -b 100 minen.pdb::A\n"
                <<argv[0]<<" c293999.xml\n";
        prf::cout<<"The option \"void_ends\" can be used to suggest the use "
                 <<"of uncharged chain ends when no capping groups are specified. \n\n";
        prf::cout<<"The option \"transformations_file\" can be used to apply "
                <<"a set of transformations to the degrees of freedom before "
                <<"calculating energy. See documentation for more information. \n";

        return 0;
    }

    int loglev=5;

    if (par.option_given("ll")) loglev=atoi(par.option("ll").c_str());
    prf::Logger::verbosity=loglev;

    string sel, filnm,extn,ifile="minen.pdb";

    if (par.n_spare_args() >=1) ifile=string(par.spare_args(0));
    prf_utils::analyze_filename(ifile,filnm,extn,sel);

    list<InstructionString> cmds=par.get_options();

    if (extn=="xml" or extn=="pdb") cmds.push_back(InstructionString("set_population "+filnm));
    else if (extn=="tconf") {
        cmds=par.get_options();
        if (ifile.substr(0,7)!=std::string("file://")) ifile="file://"+ifile;
        cmds.push_back(InstructionString("init_config "+ifile));
    } else {
        prf::cerr<<"Unknown input format for "<<argv[0]<<"\n";
        prf::cerr<<"File extension used was "<<extn<<"\n";
    }

    ph.parseCommands(cmds,argc,argv);
    if (extn=="tconf") {
        ph.init_pop();
        ph.init_coords();
        ph.reconstruct();
    }
    fh.parseCommands(cmds,argc,argv);
    fh.init_ff();
    fh.set_population(ph.population());
    ForceField *ff=fh.interaction_potential();
    ff->init();
    // Set up measurements
    Logger()(10)<<"Initializing observables and ObsHandler..\n";
    ObsHandler H;
    H.setPopulation(ph.population());
    H.set_block_props(1,"temperature");
    H.make_obs(InstructionString("RCBin HelixContent default_limits helix"));

    H.make_obs(InstructionString("RCBin BetaStrandContent default_limits strand"));H.parseCommands(cmds,argc,argv);

    // Some energy terms can depend on things calculated during some measurement
    for (size_t deptrm=0;deptrm<fh.n_obs_dependent_terms();++deptrm) {
        H.fix_dependencies(fh.obs_dependent_term(deptrm));
    }

    H.initialize();

    highlight("Energy terms for ProFASi force field "+ff->Name());
    double etot=0,eterm=0;
    for (size_t i=0;i<ff->n_terms();++i) {
        prf::cout<<ff->term(i)->Name()<<"\t\t=\t\t"
                <<(eterm=ff->term(i)->evaluate())<<"\n";
        etot+=eterm;
    }
    prf::cout<<"Total\t\t=\t\t"<<etot<<"\n";

//    ph.population()->SaveSnapshot(1,"output.pdb",0,0,etot);
    return 0;
}

/**
\page prf_energies Calculating energies of a PROFASI generated PDB/XML file
The utility \tt prf_energies calculates and prints the various PROFASI energy
terms for a PROFASI generated PDB or XML structure file. It is useful if you
have extracted a structure from a PROFASI run and want to see what the different
energy terms are like for that particular structure. When using a PDB file, you
can also use the PROFASI selection rules described in \ref prf_sel_fils to find
the energy of a subsystem (excluding everything else). \n

\section usage Usage
<tt>prf_energies structure.pdb::selections</tt>\n\n
<tt>prf_energies structure.xml</tt>\n\n

ProFASi style selection rules for PDB files may be used. For example, \n\n
to print the energy terms for a structure abc.pdb, \n
<tt>prf_energies abc.pdb</tt>\n\n
To print the energy terms involving only the atoms of chain B residues 2-10 of
abc.pdb, \n
<tt>prf_energies abc.pdb:B,2,10</tt>

To use a ProFASi text configuration file as an input, you need to set the
population in some way. This can be done, either by using the population
specifiec command line arguments or by using a settings file. This program does
not look for a settings file, unless you explicitly mention its name using the
"-st" option. In most situations, the most convenient way is to use
"--add_chain ... " through the command line.

<tt>prf_energies minen.tconf --add_chain 30 "ACE*KLVFFAE*NH2" --box_length 100</tt>\n\n

In the above example, the text configuration file minen.tconf was interpreted as
the structure information for 30 chains of sequence Acetyl--KLVFFAE--NH2. If the
sequence is long, and you would rather have it read from the original settings
file, do this:\n\n

<tt>prf_energies minen.tconf --settings_file original_settings.cnf </tt>\n\n


\note For a PDB file not generated by PROFASI, the output of this program might
be entirely incorrect. If the PDB file contains all atoms, including all hydrogen
atoms, the values will be approximately ok. Still, because of possible
differences of bond angle and bond lengths between the geometry assumed in
PROFASI and the geometry in the PDB file, there could be some bad values,
especially for the excluded volume term.

\note XML structure files store internal dihedral angles and a few rigid
body coordinates to much higher precision than is possible in a PDB file. The
energy values obtained for an XML structure file should therefore match the
ones seen during the more or less exactly, for the same force field. We use
this, to conveniently check for bugs every time we work on the energy class
implementations.


\section transformations Applying transformations to the angles
This program can optionally modify the degrees of freedom with a set of
given rules before using them. The reason for this feature is to provide a
convenient way to calculate energies for structures extracted from older
versions of ProFASi (older than 1.1).

We will illustrate with an example: the definition of the chi0 angle of most
residues is the dihedral angle N-Ca-Cb-Cg. But in the first versions of ProFASi
it used to be N-Ca-Cb-1Hb. The binary and text population configuration files
obtained from such an older version will have angles measured with a different
definition. But since this definition just creates an offset, it can be easily
corrected for. One needs to specify the (i) the type of residue where the
correction should be applied, (ii) the serial number of the side chain degree
of freedom in that residue, (iii) a constant to be added to the value of that
degree of freedom. Every residue of the same type will need just one set of
rules specifying what to do with various degrees of freedom in it.  So, one
can create a file with one line for each kind of residues for which a
correction is needed: (for example)\n\n
LYS 3 0 120 1 120 2 120\n\n
The above line says, for every lysine residue, add 120 degrees to chi0,
chi1 and chi2 before interpreting the current structure.

The file containing these transformation rules are then to be passed to
the program with the "tf" option.

*/


class torsion_transforms
{
public:
    torsion_transforms();
    torsion_transforms(std::string rulefile);
    ~torsion_transforms();
    int read_rules(std::string rulefile);
    int operator()(Ligand *lg);
private:
    std::vector<std::string> taggd;
    std::vector<std::vector<std::pair<int,double> > > rules;
};

torsion_transforms::torsion_transforms() {}

torsion_transforms::torsion_transforms(std::string rulefile) {read_rules(rulefile);}

torsion_transforms::~torsion_transforms() {}

int torsion_transforms::operator()(Ligand *lg)
{
    int ntrans=0;
    std::string resnm=lg->TLC();
    size_t res=0;

    for (;res<taggd.size();++res) if (taggd[res]==resnm) break;

    if (res!=taggd.size()) {
        for (size_t i=0;i<rules[res].size();++i) {
            double oldv=lg->get_rotDof(rules[res][i].first);
            int a0=0,a1=0;
            oldv+=rules[res][i].second;
            lg->rotDof_assign_and_reconstruct(rules[res][i].first, oldv, a0, a1);
            ++ntrans;
        }
    }

    return ntrans;
}

int torsion_transforms::read_rules(std::string rulefile)
{
    if (TestFile_r(rulefile)==0) return 0;

    std::ifstream fin(rulefile.c_str());
    std::string line;
    taggd.clear();
    rules.clear();

    while (getline(fin,line)) {
        if (line.empty()) continue;

        std::istringstream ssin(line);
        std::string resnm;
        int ncorr=0, tmpind=0;
        double magcorr=0;
        std::vector<std::pair<int, double> > tmpv;
        ssin>>resnm;
        taggd.push_back(resnm);
        ssin>>ncorr;

        for (int i=0;i<ncorr;++i) {
            ssin>>tmpind;
            ssin>>magcorr;
            magcorr=UnivConstants::pi*magcorr/180.0;
            tmpv.push_back(std::make_pair<int,double>(tmpind,magcorr));
        }

        rules.push_back(tmpv);
    }

    return 1;
}
