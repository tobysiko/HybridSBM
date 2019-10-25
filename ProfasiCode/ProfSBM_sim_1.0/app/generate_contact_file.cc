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
#include <Observables/ContactFunctions.hh>
#include <fstream>

using namespace prf;
using namespace prf_utils;
using namespace std;

int main(int argc, char *argv[])
{
    prf_utils::ProgArgs par;
    par.settings_use_type(2);
    par.option("output_file","o",1);
    par.option("cut_off","c",1,"(distance cut-off for the definition of a contact)");
    par.option("threshold","t",1,"(probability threshold to count a preserved contact)");
    par.option("aa_range","r",2,
               "(Print only contacts within a specified (C-style) residue range)");
    par.option("contact_type","ct",1,"(Proximity, CaContact, HBContact, HPContact)");
    par.option("minimum_sequence_separation", "m",1, "(to omit neighbouring residues from a contact)");
    par.option("selection","sl",1,"(selection, if you are using a PDB input)");
    par.option("log_level","ll",1);
    par.option("distances","d",1,"(print native distances to file)");
    par.option("xml","x",1,"(write native contacts to xml restraint file)");
    par.option("epsilon","e",1,"interaction strength");
    par.init_options(argc,argv);

    if (argc==1) {
        prf::highlight("Contact Map Generator");
        prf::cout<<"\n\nThis program creates a contact map file from a given structure for use "
                 <<"with ContactMap observables. The number of contacts shared with the native "
                 <<"structure by the instantaneous structure in the simulation is a measure "
                 <<"of the closeness of the two structures. For such observables, a list of "
                 <<"contacts of the native states must be prepared before a run, and its location "
                 <<"mentioned in the \"settings.cnf\" file for a given run. This program helps "
                 <<"create that precalculated list of contacts for the native structure.\n\n";
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] structure_filename\n\n";
        prf::cout<<"Options could be any number and combination of ...\n";
        par.write_available();
        prf::cout<<"\nContact_type could be any one of \n"
                 <<"{Proximity,HvContact,CaContact,HBContact,HPContact}: \n\n"
                 <<"Two residues A and B are in contact if ...\n"
                 <<"(Proximity): there exist at least two pairs of atoms (a1,b1) and (a2,b2) such that "
                 <<"a1,a2 are in A and b1,b2 are in B, with the distance between each pair "
                 <<"being less than a given cutoff. \n\n"
                 <<"(HvContact): at least one heavy atom of A is within a cutoff distance "
                 <<"of a heavy atom of B. \n\n"
                 <<"(CaContact): the C-alpha atoms of A and B are within a cutoff distance\n\n"
                 <<"(HBContact): there is a hydrogen bond between the backbone NH dipole of "
                 <<"A and backbone CO dipole of B. \n\n"
                 <<"(HPContact): the hydrophobic contact fraction between the side chains of "
                 <<"A and B is greater than a specified cutoff.\n\n";
        prf::cout<<"The structure input type is inferred from the structure file extension. "
                <<"For PDB, XML and text conf formats, use extensions \"pdb\",\"xml\" and \"tconf\"\n\n";
        prf::cout<<"If input type is pdb or xml the sequence can be inferred from the file. "
                 <<"For PDB files, you can select a range of residues to construct "
                 <<"the population from, just like in mimiqa or pdb_slices, but with"
                 <<" one difference. You have to use a separate option \"-sl\" "
                 <<"to mention the selection rather than appending it to the end"
                 <<"of the PDB file name like in mimiqa or pdb_slices. Example:\n\n"
                 <<argv[0]<<" -o contacts.dat -ct HBContact abc.pdb -sl :B,34,67\n\n"
                 <<"This is to emphasize that it uses the same selection across"
                 <<" all PDB files in the input to find the preserved contacts. "
                 <<"Example: \n\n"
                 <<argv[0]<<" -o contacts.dat -ct HBContact -sl :B,34,67 abc_*.pdb\n\n"
                 <<"Note also that the selection string is everything you would "
                 <<"append to a file name in mimiqa except for the first ':' "
                 <<"character. If you would have written file.pdb:3:A,32,48 in "
                 <<"mimiqa, or to enter the chain in a PROFASI settings file, "
                 <<"you will write 3:A,32,48. The first 3 here is the model id,"
                 <<"and can be omitted when you want to use the first model in "
                 <<"the file, as in the examples above.\n\n"

                 <<"If the input file is an XML file, at present, no selections "
                 <<"can be applied.\n\n"
                 <<"If the input is text conf file, you should provide a "
                 <<"\"settings\" file or the Population specific command line "
                 <<"instructions of BasicMCRun from which the protein sequence "
                 <<"could be determined. The name of the settings file must be "
                 <<"mentioned with the -st option. If not mentioned no settings"
                 <<" file reading is attempted.\n\n"
                 <<"Please do not confuse the option aa_range with the selection"
                 <<" system described above. Using \"-r\" does not change the "
                 <<"chain, but only omits contacts involving residues not in the"
                 <<"range. The residue range here is C-style, relative to the "
                 <<"Population in the program. The population is defined by the "
                 <<"selection system. If it is confusing, don't use the \"-r\""
                 <<"option.\n";


        return 0;
    }

    std::string sel="1:*",dummysel="";
    std::string ofile("output.nct"),ifile="somestructure.xml",xfile("output.xml");
    
    double cutv=4.5,threshold=0.75;
    std::string ctyp="CaContact",itype="pdb",filnm;
    int aa0(0),aa1(-1), loglevel=10, mssep=2;
    bool print_dist = true;
    bool contactsToXML = true;
    double epsilon = 1.0;
    double width = 1.0;
    double steepness = 1.0;
    double radius = 4.0;
    

    if (par.option_given("c")) cutv=strtod(par.option("c").c_str(),NULL);
    if (par.option_given("t")) threshold=strtod(par.option("t").c_str(),NULL);
    if (par.option_given("o")) ofile=par.option("o");
    if (par.option_given("ct")) ctyp=par.option("ct");
    if (par.option_given("m")) mssep=atoi(par.option("m").c_str());
    if (par.option_given("sl")) sel=par.option("sl");
    if (par.option_given("ll")) loglevel=atoi(par.option("ll").c_str());
    if (par.option_given("e")) epsilon=strtod(par.option("e").c_str(),NULL);
    if (par.option_given("r")) {
        aa0=atoi(par.option_arr("r",0).c_str());
        aa1=atoi(par.option_arr("r",1).c_str());

        if (aa0>aa1) std::swap(aa0,aa1);
    }
    if (par.option_given("d")) {
    	
    	if (par.option("d") == "true") print_dist = true;
    	else if (par.option("d") == "false") print_dist = false;
    	else {
    		prf::cerr<<"Invalid value for option d. Shoud be false or true.\n";
    		return 1;
    	}
    }
    if (par.option_given("x")) {
    	if (par.option("x") == "true") contactsToXML = true;
    	else if (par.option("x") == "false") contactsToXML = false;
    	else {
    		prf::cerr<<"Invalid value for option x. Shoud be false or true.\n";
    		return 1;
    	}
    }

    Logger blog;
    prf::Logger::verbosity=loglevel;
    prf::clog.open(".profasi_logfile");

    if (par.n_spare_args()==0) {
        prf::cerr<<"No structure files given.\n";
        return 1;
    } else ifile=par.spare_args(0);

    prf_utils::analyze_filename(ifile,filnm,itype,dummysel);

    PopulationHandler ph;
    list<InstructionString> cmds=par.get_options();

    if (itype=="pdb") cmds.push_back(InstructionString("set_population  "+filnm+":"+sel));
    else if (itype=="xml") cmds.push_back(InstructionString("set_population "+filnm));
    else if (itype=="tconf") {
        cmds=par.get_options();
        if (ifile.substr(0,7)!=std::string("file://")) ifile="file://"+ifile;
        cmds.push_back(InstructionString("init_config "+ifile));
    } else {
        prf::cerr<<"Unknown input format for "<<argv[0]<<"\n";
        prf::cerr<<"File extension used was "<<itype<<"\n";
    }
    
    ph.parseCommands(cmds,argc,argv);
    if (itype!="xml" and itype!="pdb") {
        ph.init_pop();
        ph.init_coords();
        ph.reconstruct();
    }
    
    
    if (!ph.population()->initialized()) {
        prf::cerr<<"generate_contact_file> Population initialization has failed."
                    <<" Can not generate contact list.\n";
        return 1;
    }

    vector<string> files;

    for (int i=0; i<par.n_spare_args(); ++i) {
        string ifile=par.spare_args(i);
        string filnm;
        prf_utils::analyze_filename(ifile,filnm,itype,dummysel);

        if (TestFile_r(filnm.c_str())) {
            files.push_back(ifile);
            blog<<"input type : "<<itype<<"\n";
        }
    }

    if (files.empty()) {
        prf::cerr<<"No valid input files.\n";
        return 1;
    }


    blog<<"output file : "<<ofile<<"\n";
    blog<<"Contact type used: "<<ctyp<<"\n";
    blog<<"Cut-off for contacts: "<<cutv<<"\n";

    Population *p=ph.population();
    std::list<SelRes> selection;

    if (aa1<0) aa1=p->NumberOfResidues();

    blog<<"Ligand index range : "<<aa0<<"(inclusive) to "<<aa1<<"(exclusive)\n";

    Proximity prx;
    CaContact cacn;
    HBContact nhb;
    HPContact hpc;
    bool assymetric_contact=false;
    ContactFunction *cf;

    if (ctyp==string("Proximity")) {
        cf=&prx;
    } else if (ctyp==string("HvContact")) {
        cf=&prx;
        prx.minimum_links(1);
    } else if (ctyp==string("CaContact")) {
        cf=&cacn;
    } else if (ctyp==string("HPContact")) {
        cf=&hpc;
    } else if (ctyp==string("HBContact")) {
        assymetric_contact=true;
        cf=&nhb;
    } else {
        std::cerr<<"Unknown contact type.\n";
        return 1;
    }
    
    
    cf->min_seq_sep(mssep);

    if (par.option_given("c")) cf->set_cutoff(cutv);

    cf->init(p);
    deque<int> i1,i2;
    
    Matrix<double> prob_intct(aa1-aa0,aa1-aa0);
    
    
    
    vector<vector < list <double> > > distances(aa1-aa0, vector< list<double> >(aa1-aa0, list<double>(0))); // NEW! -TS
    
    
    prob_intct*=0;
    list<AtomRecord> rcd;
    int nstrucs=0;
    double dist; // NEW! -TS
    
    for (size_t k=0; k<files.size(); ++k) {
        if (itype==string("pdb")) {
            PDBReader pdb(files[k]);
            
            if (pdb.read_matrix()==0) {
            	
                prf::cerr<<"Unable to read data from file "<<files[k]<<"\n";
                continue;
            }
            
            selection.clear();
            
            pdb.mk_selection(sel,selection);
            
            pdb.records(selection,rcd);
            std::vector<bool> specified(p->NumberOfAtoms(),false);
            p->ImportStructure(rcd,specified);
            
            p->guess_missing_coordinates(specified);
            
        } else if (itype=="xml") {
            ph.read_xml_pop(files[k]);
        } else {
            FILE *fp=fopen(files[k].c_str(),"r");
            p->ReadConf_text(fp);
            fclose(fp);
            ph.reconstruct();
        }
        for (int i=aa0; i<aa1; ++i) {
            for (int j=aa0; j<aa1; ++j) {
            	if (ctyp==string("CaContact")&&print_dist){  // NEW! -TS
            		dist = (*cf)(i,j,true);
            		//prf::cout<<i<<" "<<j<<" "<<dist<<"\n";
            		if (dist != -1.0) {
            			prob_intct[i-aa0][j-aa0]+=1;
            			distances[i-aa0][j-aa0].push_back( sqrt( dist ) );
            			prf::cout<<i<<" "<<j<<" "<<distances[i-aa0][j-aa0].back()<<"\n";
            		}
            	}
            	if ((*cf)(i,j)) prob_intct[i-aa0][j-aa0]+=1;
            	
            }
        }

        ++nstrucs;
    }
    
    
    if (nstrucs==0) {
        prf::cerr<<"No valid structures.\n";
        return 1;
    }

    for (int i=0; i<(aa1-aa0); ++i) {
        for (int j=0; j<(aa1-aa0); ++j) {
            prob_intct[i][j]/=nstrucs;
//            prf::cout<<(int) prob_intct[i][j]<<" ";
        }
//        prf::cout<<"\n";
    }

    for (int i=0; i<(aa1-aa0); ++i) {
        for (int j=0; (assymetric_contact&&j<(aa1-aa0))||j<i; ++j) {
            if (prob_intct[i][j]>threshold) {
                bool usectc=false;

                if (assymetric_contact || prob_intct[j][i]>threshold) usectc=true;

                if (usectc) {
                    i1.push_back(i+aa0); i2.push_back(j+aa0);
                } else
                    prf::cerr<<"Assymetry: ("<<i+aa0<<", "<<j+aa0<<") "
                            <<"were found to be a contact. But not "
                            <<"("<<j+aa0<<", "<<i+aa0<<")\n"
                             <<"This pair will be ignored. \n\n";
            }
        }
    }

    prf::cout<<"Found "<<i1.size()<<" preserved contacts in the input file(s).\n";
    
    int chain1,chain2;
    std::string resname1,resname2;
    int res1,res2;
    
    ofstream xmloutLJ("LJcontacts.xml");
    ofstream xmloutGAUSS("GAUSScontacts.xml");
    ofstream bondout("pseudobonds.txt");
    if (contactsToXML){
    
    xmloutLJ<<"<distance_restraints>\n  <formatted_data>"
    <<"\n    <format name=\"restraint\" type=\"$3\">\n      <atom1>$1</atom1>"
    <<"\n      <atom2>$2</atom2>\n    <parameters>\n      <epsilon>$4</epsilon>"
    <<"\n      <minimum>$5</minimum>\n    </parameters>\n    </format>\n    <data>\n";
    
    xmloutGAUSS<<"<distance_restraints>\n  <formatted_data>"
        <<"\n    <format name=\"restraint\" type=\"$3\">\n      <atom1>$1</atom1>"
        <<"\n      <atom2>$2</atom2>\n    <parameters>\n      <radius>$4</radius>"
        <<"\n      <steepness>$5</steepness>\n      <minimum>$6</minimum>\n      <width>$7</width>"
        <<"\n      <depth>$8</depth>\n    </parameters>\n    </format>\n    <data>\n";
    
    }
    if (i1.size()>0) {
        ofstream fout(ofile.c_str());

        for (int i=0; i<p->NumberOfChains(); ++i) {
            fout<<"# Chain "<<i<<": "<<p->PepName(i)<<"\n";
        }
        
        double d;
        for (size_t i=0; i<i1.size(); ++i) {
        	d = distances[ i1[i] ][ i2[i] ].back();
        	res1 = i1[i];
	        res2 = i2[i];
	        chain1 = 0;
	        chain2 = 0;
	        AminoAcid* am1 = p->Chain(chain1)->AA(res1);
	        AminoAcid* am2 = p->Chain(chain2)->AA(res2);
	        resname1 = am1->TLC();
	        resname2 = am2->TLC();
	        
	        Atom* atm1=p->ligand(res1)->labeled_atom(" CA ");
	        Atom* atm2=p->ligand(res2)->labeled_atom(" CA ");
	        
	        if (contactsToXML){
	        
	        
	        	
	        	xmloutLJ << "        "<<chain1<<"/"<<res1<<"/" << resname1 << "/_CA_ " <<chain2<<"/"<<res2<<"/"<< resname2 << "/_CA_ " << " LJ "<< epsilon << " " << d << "\n";
	        	xmloutGAUSS << "        "<<chain1<<"/"<<res1<<"/" << resname1 << "/_CA_ " <<chain2<<"/"<<res2<<"/"<< resname2 << "/_CA_ " << " GaussContact "  << radius << " "<< steepness << " " <<  d <<" "<< width <<" "<<epsilon<< "\n";
	        	// 0/0/GLY/_CA_ 0/15/GLU/_CA_ quadratic 5.5 3.0
	        	// 0/1/GLU/_CA_ 0/14/THR/_CA_ quadratic 5.5 3.0
	        	// 0/2/TRP/_CA_ 0/13/VAL/_CA_ quadratic 5.5 3.0
	        	// 0/3/THR/_CA_ 0/12/THR/_CA_ quadratic 5.5 3.0
	        	// 0/4/TYR/_CA_ 0/11/PHE/_CA_ quadratic 5.5 3.0
	        	prf::cout << atm1->UniqueId() << " " << atm2->UniqueId() << " "<< res1 << " " << res2 << "\n";
	        	bondout << "#0:" << res1 << "@ca " << "#0:" << res2 << "@ca blue "<< d << "\n";
	        }
	        
	        if (print_dist){
	        	
	        		fout<<i<<"  "<<i1[i]<<"  "<<i2[i]<<"  "<< d <<"\n";
	        	
	        }
	        else{
	        	
	        		fout<<i<<"  "<<i1[i]<<"  "<<i2[i]<<"\n";
	        	
	        }
        }

        fout.close();
        
    }
    
    if (contactsToXML){
        
    xmloutLJ<<"    </data>\n  </formatted_data>\n</distance_restraints>\n";
    xmloutGAUSS<<"    </data>\n  </formatted_data>\n</distance_restraints>\n";
    }
    
    xmloutLJ.close();
    xmloutGAUSS.close();
    bondout.close();
    prf::cout<<"Contact type "<<ctyp<<", ";
    prf::cout<<"cutoff = "<<cf->get_cutoff();

    if (ctyp==string("Proximity") || ctyp==string("HvContact") ||
        ctyp==string("CaContact")) {
        prf::cout<<" Angstroms.\n";
    } else if (ctyp==string("HBContact")) {
        prf::cout<<" PROFASI energy units. (Note: HBContact is asymmetric!)\n";
    } else if (ctyp==string("HPContact")) {
        prf::cout<<" of maximum hydrophobic contact between residue pairs.\n";
    }

    return 0;
}

/**
  \page contactgen Generating native contact lists
  The program \c generate_contact_file can be used to generate a list of
  contacts found in a structure. Such a file can be used to initialize a
  prf::ContactMap observable, to monitor the on/off state of native contacts
  through a run.

  \section contactgen_1 Types of contact maps in PROFASI
  For the purpose of contact map measurements, "contact" in PROFASI is between
  two residues. Therefore it can be enumerated by two integers corresponding to
  their unique serial numbers in the population. The following types of contacts
  are available by default:
  \li \e CaContact : the \f$C_\alpha\f$ atoms of the residues are within a cut-off
  distance
  \li \e HvContact : at least one heavy atom in one residue is within a cut-off
  distance of a heavy atom in the other residue
  \li \e Proximity : there exist at least two pairs of heavy atoms, such that
  atoms in a pair do not come from the same residue, and are within a cut-off
  distance
  \li \e HBContact : there is a backbone hydrogen bond from the NH of the left-hand
  residue to the CO of the right-hand residue. Note that it is asymmetric!
  \li \e HPContact : there is a hydrophobic contact between the residues

  \section contactgen_2 Typical usage
  Given a PDB file \c abc.pdb, you can find the hydrogen bond contacts in it
  like this:

  \verbatim
  generate_contact_file -ct HBContact -o result.dat abc.pdb
  \endverbatim

  If the PDB file contains the structure of a large protein and you are
  interested in the internal contacts of a segment between residues 50--100,
  you can do this:

  \verbatim
  generate_contact_file -ct HBContact -o result.dat abc.pdb -sl :A,50,100
  \endverbatim

  In the above, we gave the selections differently compared to PROFASI programs
  like mimiqa. To select residues 50 to 100 in chain A of file \c abc.pdb in
  \c mimiqa, we would write \c abc.pdb:1:A,50,100. The syntax for the selections
  for \c generate_contact_file is in fact the same. But one has to specify the
  selection string separately with the \c -sl option. Think of the mimiqa
  filename-selection combination as \c filename:selection_string. For
  \c generate_contact_file, we need the \c selection_string part. That would
  have been \c 1:A,50,100 above, but if we are happy to use the first model
  in the PDB file, we can omit the model identifier \c 1 here. This explains
  the string \c :A,50,100 above.

  The reason for the separate option \c -sl for selections is that it can be
  applied to a bunch of PDB files together. The program then looks for preserved
  contacts in all those files. Something like this:
  \verbatim
  generate_contact_file -ct HBContact -o result.dat abc_*.pdb -sl :A,50,100
  \endverbatim

  \note The directionality of the hydrogen bond potential in PROFASI means
  that the HBContact measure may miss many hydrogen bonds in a PDB file,
  typically refined with some other force field. To find all hydrogen bonds, it
  is better to first regularize the structure (See \ref regul), and take the
  resulting energy minimized structure file as input. The XML output files of
  the regularizer may be converted to a PDB file with prf_convert if you want
  to use selections like in the examples above.

  \section contactgen_exls Example with output
  Let's take a concrete example: the C-terminal hairpin of protein G, and
  generate the native hydrogen bond contacts. This hairpin has been studied
  as an excised peptide in experiments and numerous simulations. Download
  1GB1.pdb from the PDB. The hairpin consists of residues 41--56 of the only
  chain in this structure.

  First, let's cut out the relevant residues into a new PDB file.
  \verbatim
  pdb_slices 1GB1.pdb::A,41,56 -o hairpin.pdb
  \endverbatim
  Check the contents of the file hairpin.pdb. It should only have the atom
  records for residues 41 through 56. Next, let's regularize this file:
  \verbatim
  regularize hairpin.pdb
  \endverbatim
  There should be two output files, "min_etot.xml" and "min_rmsd.xml". We need
  the first, to get the native hydrogen bonds.
  \verbatim
  generate_contact_file -ct HBContact min_etot.xml -o hairpin.nhb
  \endverbatim
  This will generate an output file \c hairpin.nhb with the following contents:
  \verbatim
# Chain 0:   * GEWTYDDATKTFTVTE *
0  1  14
1  3  12
2  5  10
3  8  5
4  9  5
5  10  5
6  12  3
7  14  1
  \endverbatim
  The first line shows the sequence of the peptide segment. The residue indexes
  in the bond listing must be interpreted relative to this sequence, and they
  start from 0. So, the hydrogen bonds identified were 2 each between residue
  pairs 1(E)-14(T), 3(T)-12(T), 5(D)-10(T) and another between the NH dipole
  of 8(T) and the CO dipole of 5(D). PROFASI's regularizer uses a little bit of
  MC. So, you could run the regularizer a few times and choose the lowest energy
  regularized state found.

  \section contactgen_opts Options supported by generate_contact_file

    \li \b --output_file or \b -o : Name of the output file
    \li \b --contact_type or \b -ct : Contact type. Possible values: Proximity,
CaContact, HBContact, HPContact
    \li \b --selection or \b -sl : Select a certain range in a PDB file and work
    on it.
    \li \b --aa_range or \b -r : Exclude contacts outside a specified (C-style)
residue range
    \li \b --minimum_sequence_separation or \b -m : Exclude pairs within a
certain sequence separation
    \li \b --cut_off or \b -c : Distance cut-off for the definition of a contact
    \li \b --threshold or \b -t : Probability threshold to count a preserved
 contact
    \li \b -log_level or \b -ll : Log level of output messages

    Some of these options have been described above in the examples, and the
    others have fairly obvious meanings. But we would like to point out that
    the option \c --aa_range or \c -r is \b not the same as the option
    \c --selection . When you use \c -sl, you indicate that you want to work
    with the selected part of the protein in the simulation. When you use
    \c -r you indicate that you want a measurement of a limited part of all
    native contacts of the simulated system.

    For example, let's consider the situation when we want to simulate the whole
    protein in 1GB1.pdb. We want to monitor only the hydrogen bonds corresponding
    to the helix (23 (ALA) to 36 (ASP)). The required contacts file will then
    be generated by
    \verbatim
    generate_contact_file -ct HBContact 1GB1_reg.pdb -o 1GB1_hel.nhb -r 22 35
    \endverbatim

    In the above, 1GB1_reg.pdb is obtained by first regularizing 1GB1.pdb and
    then converting the resulting min_etot.xml to a PDB file. The file
    1GB1_hel.nhb will contain only a subset of all native contacts in 1GB1.pdb,
    but will number them relative to the sequence of whole chain. If you used
    the \c -sl option, the same contacts will be generated, but they will be
    numbered relative to the segment of the helix section alone. If that's what
    you want to simulate, fine. But in the present example, those indexes will
    correspond to the wrong residues of the chain.

    \note The residue ranges in the \c -r option start from 0, while the \c -sl
    option makes a selection on the PDB file with numbering as in that file.

    \sa prf::ContactMap, \ref regul, \ref prf_convert, \ref pdb_slices

  */
