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

#include <Algorithms/Minimizer.hh>
#include <Aux/ProgUtils.hh>
#include <Aux/profasi_version.hh>
#include <Aux/prf_time.hh>
#include <Elements/PopulationHandler.hh>
#include <Energy/FFHandler.hh>
#include <Observables/ProteinRMSD.hh>
#include <Observables/ObsExpression.hh>
#include <Aux/RandomNumberHandler.hh>
#include <Algorithms/MC.hh>
#include <Updates/Rot.hh>
#include <Updates/BGS.hh>
#include <sstream>

using namespace prf;
using namespace prf_utils;

/**
  \page regul Regularization: approximating a protein structure
  \section regul_def Definition
  Regularization is the process of identifying the best approximation of a given
  protein structure which, (i) satisfies the constraints of the protein model,
  such as the bond length and bond angle values imposed by PROFASI (ii) is a
  minimum of the interaction potential.
  \section regul_expl Explanation
  Although the protein model in PROFASI is an all atom model, it works under
  the approximation that the bond lengths and bond angles do not change. Further,
  these geometrical properties are parameters of the model. The parameters of
  the energy function are derived assuming certain values for the bond lengths
  and bond angles. For instance, the \f$N-C_{\alpha}\f$ bond for any residue in
  PROFASI is 1.46 &Aring;. This an approximation derived from a statistical
  analysis of structures in the PDB.

  A typical protein structure downloaded from the PDB will in general not have
  a \f$N-C_{\alpha}\f$ bond of exactly 1.46 &Aring;, but something close to
  it. Using the atom coordinates as given in the PDB file will give incorrect
  energies in PROFASI, because the values of the bond lengths and angles are
  implicit in the derivation of parameters of our energy function. We therefore
  need to find alternative coordinates for the atoms, as close as possible to
  the ones in the PDB structure, which satisfy the bond length and bond angle
  constraints. This can be done by minimizing the all-atom RMSD with the given
  protein structure with respect to the degrees of freedom in PROFASI.

  But the structure obtained from the RMSD minimization may be of limited value.
  The structure files in the protein data bank are normally refined using some
  force fields. The precise location of the atoms which leads to the lowest
  energy in the refining force field will not in general lead to a reasonably
  low energy in PROFASI's force field. If simulations are to be started using
  an approximated structure, or the energy calculated is to be compared to
  energies obtained in simulation, it is necessary to find a structure close to
  the given structure which is a minimum of PROFASI's force field. A direct
  minimization of energy after reading in the coordinates from a PDB file
  generally leads to rapid breaking of a lot of secondary structure. This is
  because using torsion angles calculated from the Cartessian coordinates in
  the PDB files to initialize a protein structure often leads to a clashes
  between atoms somewhere in the structure. Clashes have such high energies
  in PROFASI that it is then seen as a reduction in energy if we move the
  clashing atoms apart even if we break secondary structure in the process. For
  this reason, the energy minimization is done in many stagees with a gradually
  decreasing RMSD constraint. The program \e regularize performs this using
  a mixture of Monte Carlo and conjugate gradient minimization.

  \section regul_howto How to regularize structures in PROFASI
  Given a PDB file abc.pdb, to find a regularized structure, do the following:
  \verbatim
  regularize abc.pdb
  \endverbatim
  This will try to find a regularized structure using default values for a set
  of options. To see the various supported options controlling the behaviour
  of the program, see the documentation of \ref regularizer .

  Two output files are generated. \e min_rmsd.xml is the best
  approximation found by minimizing RMSD alone. \e min_etot.xml corresponds to
  a minimum of PROFASI's energy function near the given structure. The output
  files are in \ref xmlstruc. They can be converted to the PDB
  format like this:
  \verbatim
  prf_convert min_etot.xml min_etot.pdb
  \endverbatim

  \section regul_caveats Caveats
  Both the processes of minimizing RMSD and energy might fail to produce good
  approximations for different reasons. When RMSD minimization fails to get
  a satisfactorily low value, it is most often because in those cases atoms
  are labeled differently in PROFASI and in the PDB file. Labeling errors in
  PDB files are unfortunately not infrequent. This results in the program
  trying to bring the wrong atoms close to each other. Less frequently, there
  are structures in which an amino acid at a given position might exhibit bond
  lengths or bond angles which are relatively far from the approximations (based
  on average values over PDB) used in PROFASI. The minimization process then
  finds a compromise which may not look like a good approximation.

  Since the "landscape" of the interaction potential as a function of the
  degrees of freedom is much rougher than the lanscape of RMSD, it is much more
  likely that the energy minimization process gets trapped in an uninteresting
  local minimum of high value. The minimum of energy closest to the RMSD
  minimized structure is not necessarily interesting, as it may still contain
  clashes. The lowest possible energy within 2 &Aring; of the RMSD minimum may
  be the 23427th minimum ranked by RMSD. The "best" structure is therefore a
  nebulous concept. The program \e regularize therefore has a more modest goal
  of producing "a" minimum of energy, rather than "the" minimum. Quite often
  the minimum produced is also of low RMSD and energy. Since it uses
  Monte Carlo as part of its minimization algorithm, one should run it a few
  times and select one of the output structures. If the protein folds with
  PROFASI's force field, the full length folding simulation will almost
  certainly find structures with lower energy than what \e regularize finds.

  We should also mention that sometimes, regularization does not preserve all
  parts of the given structure. It may be that the structure under consideration
  is not even a local minimum in the energy function of PROFASI, so that the
  energy minimization process will move away from it. One has to remember that
  the energy function is a work in progress.

  \sa regularizer
  */

//! A program to regularize PDB structures
/**
  The goals of the program \e regularize are explained in detail in \ref regul.
  Usage:
  \verbatim
  regularize [OPTIONS] somefile.pdb
  \endverbatim
  where options could be :
  \li \b --geometric_only or \b -g : Stop after RMSD minimization
  \li \b --good_enough_rmsd or \b -gr : Stop conjugate gradient minimization
  cycles for the RMSD minimization part once RMSD falls below a given value.
  \li \b --num_e_minim_cycles or \b -nec : Number of cycles of minimization
  of E+k*RMSD . A cycle with k=0 is always done at the end.
  \li \b --initial_RMSD_factor or \b -f0 : Initial factor k in e+k*RMSD
  \li \b --final_RMSD_factor or \b -f1 : Final factor k in e+k*RMSD
  \li \b --log_level or \b -ll : Verbosity of log messages
  \li \b --mccyc_per_minim or \b -nmc : Number of MC cycles before each
  conjugate gradient minimization cycle
  \section prog_reg_examples Examples
  The simplest use case is
  \verbatim
  regularize abc.pdb
  \endverbatim

  This produces two files \e min_etot.xml and \e min_rmsd.xml.

  \verbatim
  regularize -g -gr 0.5
  \endverbatim
  The above will produce only a file \e min_rmsd.xml and stop either when the
  RMSD minimization process converges or when the RMSD value reaches 0.5 &Aring;.

  \sa \ref regul
  */
class regularizer {
public:
    regularizer();
    ~regularizer();
    int init(int argc, char *argv[]);
    int execute();
    int minimize_rmsd();
private:
    ProgArgs prog;
    PopulationHandler PH;
    Population *p;
    ProteinRMSD rmsd;
    FFHandler ffh;
    Minimizer cg;
    ObsExpression m;
    std::string rminfile, erminfile, ifile;
    int mccycperminim,nstages;
    double rmsfact0,rmsfact1;
};

regularizer::regularizer()
{
    prog.settings_use_type(2);
    prog.option("log_level","ll",1,"(verbosity level of log messages)");
    prog.new_switch("geometric_only","g","false","(Minimize RMSD and stop)");
    prog.option("initial_structure","start",1,
                "(where to start, if different from the one being regularized)");
    prog.option("num_e_minim_cycles","nec",1,
                "(number of minimization cycles of E+k*RMSD)");
    prog.option("initial_RMSD_factor","f0",1,
                "(initial factor k in E+k*RMSD)");
    prog.option("final_RMSD_factor","f1",1,
                "(final factor k in E+k*RMSD)");
    prog.option("good_enough_rmsd","gr",1,
                "(abort RMSD minimization below a given value)");
    prog.option("mccyc_per_minim","nmc",1,
                "(number of MC cycles before each CG minimization)");
    p=NULL;
    mccycperminim=10;
    nstages=9;
    rmsfact0=256.0;
    rmsfact1=1.0;
}

regularizer::~regularizer() {}

int regularizer::init(int argc, char *argv[])

{
    prog.init_options(argc,argv);

    if (prog.n_spare_args()==0) {
        prf::cout<<"This program takes a PDB file as input and produces \n"
        <<"two output files. The first represents the best approximation \n"
        <<"of the given PDB file that can be generated within the bond \n"
        <<"length and bond angle constraints of ProFASi, using conjugate \n"
        <<"gradient minimization of RMSD. The second output \n"
        <<"structure is represents an energy minimum of the force field \n"
        <<"in the neighbourhood of the starting structure. Both output \n"
        <<"files are written in ProFASi's XML structure format.\n"
        <<"Usage: "<<argv[0]<<" [OPTIONS] something.pdb\n"
        <<"Where OPTIONS could be one or more of ...\n";
        prog.write_available();
        return 0;
    }

    ifile=prog.spare_args(0);
    rminfile="min_rmsd.xml";
    erminfile="min_etot.xml";

    int loglevel=20;
    if (prog.option_given("ll")) loglevel=atoi(prog.option("ll").c_str());
    prf::Logger::verbosity=loglevel;

    if (prog.option_given("gr")) {
        double rmsdthr=strtod(prog.option("gr").c_str(),NULL);
        cg.set_good_enough(rmsdthr);
    }

    if (prog.option_given("nec")) {
        nstages=atoi(prog.option("nec").c_str());
    }

    if (prog.option_given("nmc")) {
        mccycperminim=atoi(prog.option("nmc").c_str());
    }

    if (prog.option_given("f0")) {
        rmsfact0=strtod(prog.option("f0").c_str(),NULL);
    }
    if (prog.option_given("f1")) {
        rmsfact1=strtod(prog.option("f1").c_str(),NULL);
    }
    std::list<InstructionString> cmds;

    if (prog.option_given("start")) {
        std::string stfile=prog.option("start");
        cmds.push_back(InstructionString("set_population "+stfile));
    } else cmds.push_back(InstructionString("set_population "+ifile));

    PH.parseCommands(cmds,argc,argv);
    p=PH.population();
    if (PH.init_pop()==0) {
        prf::cerr<<"Population initialilsation failed. Population can be set up "
                <<"using command line arguments or the settings file. \n";
        prf::cerr<<"\nRun "<<argv[0]<<" -h \n\nfor usage information.\n";
        return 0;
    }
    PH.init_coords();
    PH.reconstruct();
    m.set_population(p);

    ffh.parseCommands(cmds,argc,argv);

    return 1;
}

int regularizer::minimize_rmsd()
{
    Logger blog(10);
    rmsd.setPopulation(p);
    rmsd.set_logger_threshold(30);
    rmsd.set_struc1(ifile);
    rmsd.live_first_struc(true);
    rmsd.filter("+all");

    for (int i=0;i<p->num_chains();++i) {
        std::ostringstream ost;
        ost<<ifile<<"::#"<<i;
        rmsd.set_struc1(ost.str());
        ost.str("");
        ost<<"$::"<<(char)('A'+i);
        rmsd.set_struc2(ost.str());
        ost.str("");
        ost<<"FullRMSD_"<<(char)('A'+i);
        rmsd.Name(ost.str());
        rmsd.init();
        m.reset_dof_list();

        for (int idof=0;idof<p->Chain(i)->n_dof();++idof) {
            m.add_dof(p->get_dof_info(i,idof));
        }
        m.reset_obs_list();
        m.add_obs(1.0,&rmsd);
        m.init();
        m.value();

        cg.set_function(&m);
        prf::cout<<"Minimizing "<<rmsd.Name()<<"\n";
        cg.minimize();
    }

    rmsd.set_struc1(ifile);
    rmsd.set_struc2("$::*");
    rmsd.Name("FullRMSD_all_chains");
    rmsd.init();
    m.reset_dof_list();

    for (int idof=0;idof<p->n_dof();++idof) {
        m.add_dof(p->get_dof_info(idof));
    }
    m.reset_obs_list();
    m.add_obs(1.0,&rmsd);
    m.init();

    cg.set_function(&m);
    cg.minimize();

    blog<<"Saving RMSD minimized file...\n";
    FILE *fp=fopen(rminfile.c_str(),"w");
    fprintf(fp,"<?xml version=\"1.0\"?>\n<structure>\n");
    fprintf(fp,"<profasi_version>%s</profasi_version>\n",
            profasi_version().c_str());
    fprintf(fp,"<creation_time>\nUTC %s</creation_time>\n",
            prf_time().to_UTC().c_str());
    fprintf(fp,"<box_length>%.16f</box_length>\n",AtomCoordinates::boxL());
    fprintf(fp,"<remark>\n");
    fprintf(fp,"Structure with the smallest all atom RMSD with ");
    fprintf(fp,"%s within ProFASi's bond length and bond angle ",ifile.c_str());
    fprintf(fp,"constraints.\n</remark>");
    p->Write_XML(fp);
    fprintf(fp,"</structure>\n");
    fclose(fp);
    return 1;
}

int regularizer::execute()
{
    Logger blog(10);

    minimize_rmsd();
    if (prog.switch_given("g")) return 0;

    blog<<"Initialising force field...\n";
    ffh.init_ff();
    prf::ObsEnergy oen;
    oen.setObsName(rmsd.Name());
    oen.connectObs(&rmsd);
    ffh.useEnergy(&oen);
    ffh.set_population(p);
    ffh.interaction_potential()->init();
    double minen;
    blog<<"Initial value of energy = "<<
            (minen=ffh.interaction_potential()->evaluate())<<"\n";

    blog<<"Starting energy minimization...\n";

    RandomNumberHandler rh;
    rh.auto_seed(0);
    prf::MC mc;
    mc.Connect(PH.population());
    mc.forcefield_handler(&ffh);
    mc.RandomNumberGenerator(rh.generator());
    mc.updates_hander()->set_n_temps(1);

    mc.Setup();
    mc.temperature_in_kelvin(274);

    double rmsstep=0.5,rmsfactor=rmsfact0;
    if (nstages!=1) {
        rmsstep=std::pow(rmsfact1/rmsfact0,1.0/(nstages-1));
    }


    cg.set_abs_scale();

    for (int i=0;i<nstages+1;++i) {
        if (i==nstages) rmsfactor=0;
        prf::cout<<"Stage "<<i<<"\n";
        prf::cout<<"Initial value of pseudo-energy = "
                <<ffh.interaction_potential()->value()<<"\n";
        prf::cout<<"Reducing RMSD factor to "<<rmsfactor<<"\n";
        oen.setScaleFactor(rmsfactor);
        double emin=ffh.interaction_potential()->evaluate();
        prf::cout<<"New value of pseudo-energy after reduction = "<<emin<<"\n";
        std::vector<double> mincoords;
        p->get_dof(mincoords);
        prf::cout<<"Running "<<mccycperminim<<" MC cycles looking for a lower "
                <<"energy than the starting value "<<emin<<"...\n";
        for (int icyc=0;icyc<mccycperminim;++icyc) {
            mc.RunCycle();
            if (ffh.interaction_potential()->value()<emin) {
                p->get_dof(mincoords);
                emin=ffh.interaction_potential()->value();
                prf::cout<<"New lowest function energy "<<emin<<"\n";
            }
        }
        prf::cout<<"Restoring system to the best known conformation for "
                <<rmsfactor<<"*RMSD+physical_energy = "<<emin<<"\n";
        p->set_dof(mincoords);
        p->Reconstruct();
        ffh.interaction_potential()->evaluate();

        prf::cout<<"Starting CG minimisation loops...\n";
        m.reset_obs_list();
        for (size_t j=0;j<ffh.interaction_potential()->n_terms();++j) {
            m.add_obs(1.0,ffh.interaction_potential()->term(j));
        }
        m.init();

        cg.minimize();
        ffh.interaction_potential()->refresh();
        double lasten=ffh.interaction_potential()->value();
        if (lasten<minen) {
            prf::cout<<"New lowest energy structure. \n";
            ffh.interaction_potential()->print_contributions(prf::cout);

            blog<<"Saving minimum energy file...\n";
            FILE * fp=fopen(erminfile.c_str(),"w");
            fprintf(fp,"<?xml version=\"1.0\"?>\n<structure>\n");
            fprintf(fp,"<profasi_version>%s</profasi_version>\n",
                    profasi_version().c_str());
            fprintf(fp,"<creation_time>\nUTC %s</creation_time>\n",
                    prf_time().to_UTC().c_str());
            fprintf(fp,"<box_length>%.16f</box_length>\n",AtomCoordinates::boxL());
            fprintf(fp,"<remark>\n");
            fprintf(fp,"An energy minimum structure with ProFASi's force field ");
            fprintf(fp,"close to the given input file %s\n",ifile.c_str());
            fprintf(fp,"</remark>");
            p->Write_XML(fp);
            fprintf(fp,"</structure>\n");
            fclose(fp);
            minen=lasten;
        }
        rmsfactor*=rmsstep;;
    }

    return 0;
}

int main(int argc, char *argv[])
{
    regularizer reg;
    if (reg.init(argc,argv)) {
        return reg.execute();
    }
    return 0;
}
