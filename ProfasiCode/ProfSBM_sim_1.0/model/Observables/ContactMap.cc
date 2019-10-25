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

#include "ContactMap.hh"
#include <fstream>
#include <sstream>
#include "../Aux/profasi_io.hh"

using std::bitset;
using std::string;
using std::ifstream;
using std::istringstream;
using std::vector;

using namespace prf;

using namespace prf_utils;

ContactMap::ContactMap()
{
    n_monitored_contacts=0;representation_size=1;
    mainmap.reset();cmap.resize(representation_size);
    contact.clear();
    Name("ContactMap");
    mycf=false;
    fixed_his=true;
    contact_function=NULL;
}

ContactMap::~ContactMap() {}

void ContactMap::his_range(double xmn, double xmx) {}

void ContactMap::his_nbins(int n) {}

void ContactMap::rangeEstimate(double &x0, double &x1)
{
    x0=-0.5;x1=n_monitored_contacts+0.5;
    if (!userbinsz) xbin0=1;
}

/**
\page opt_ContactMap ContactMap
\section options Available options
<ul>
<li><b>type</b> type contact_type <br>
contact_type can be CaContact, for C-alpha contacts, HvContact, for heavy-atom
contacts, HBContact, for hydrogen bond contacts, Proximity, for a two-link
 closeness measure, HPContact for hydrophobic contacts.
</li>
<li><b>file</b> file contactmapfile<br>The file contactmapfile should contain 3
 columsns and list all contacts to be monitored. The first column is just an
 enumerator. The second and third columns should contain integers identifying the
 contact. A contact is always a contact between two objects, and the two integers
 must identify those two objects to the contact type. In all types of
 ContactFunctions in PROFASI, the integers represent the Ligand_indices of the
 residues relative to the entire system. </li>
<li><b>opt</b> opt cutoff_values <br>There could be multiple cut-off values that
 make sense for a contact function. So, a generic "opt" option is given so that
 more than one such value can be specified one after the other. For CaContact,
 there is only one cutoff: the maximum separation between two C-alpha atoms for
 them to be considered to be in contact. For HBContact, there is the minimum
 hydrogen bond strength between two residues. For HPContact, there is a
 hydrophobic contact fraction upon which a cut-off should be used. For heavy atom
 contacts and for Proximity, distance. For any kind of ContactFunction, one can
 put an additional cut-off on the minimum sequence separation for a pair to have
 a contact.</li>
\section examples Examples
new_obs ContactMap NHB type HBContact ; file native_hb.list<br>
This creates a new ContactMap observable to track hydrogen bond contacts between
 the pairs listed in file "native_hb.list". The measurements are going to be
 referred to as "NHB" in the run.  <br><br>
new_obs ContactMap CaC type CaContact ; file native_cac.list ; opt cutoff 6.0<br>
This creates a new ContactMap observable to track C-alpha contacts between the
 pairs listed in file "native_cac.list", using a cutoff of 6.0 Angstroems. The
 "opt cutoff 6.0" can be omitted. <br><br>
new_obs ContactMap HPC type HPContact ; file native_hp.list <br>
A hydrophobic contact monitor.
</ul>

\sa prf::ContactMap
*/
int ContactMap::init_obs()
{
    if (Observable::init_obs()==0) return 0;

    std::vector<std::string> parts,contopts;

    std::string cftype="unknown";

    //The contact type has to be set explicitly by a user command
    //The rest of the processing depends on it. So, the commands determining
    //and configuring the contact function will be processed first.
    for (size_t i=0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);
        if (parts[0]==string("type") && parts.size()>=2) {
            string tp=parts[1];

            if (tp=="CaContact" or tp=="Proximity" or tp=="HvContact"
                or tp=="HBContact" or tp=="HPContact") cftype=tp;
            else {
                prf::cerr<<"Unknown contact type "<<tp<<" requested.\n";
                prf::cerr<<"Recognized contact types are \n"
                <<"CaContact : for C-alpha contacts\n"
                <<"HvContact : for heavy atom contacts\n"
                <<"Proximity : for proximity contacts\n"
                <<"HBContact : for hydrogen bond contacts\n"
                <<"HPContact : for hydrophobic contacts\n";
            }
        }
        if (parts[0]=="opt" && parts.size()>=2) {
            string copt="";
            for (size_t j=1;j<parts.size();++j) copt+=(parts[j]+string(" "));
            contopts.push_back(copt);
        }
    }

    if (cftype==string("unknown")) return 0;

    mycf=true;

    if (cftype==string("CaContact")) contact_function=new CaContact();

    if (cftype==string("HvContact") or cftype==string("Proximity")) {
        Proximity *tmp=new Proximity();

        if (cftype=="HvContact") tmp->minimum_links(1);

        contact_function=tmp;
    }

    if (cftype==string("HBContact")) contact_function=new HBContact();

    if (cftype==string("HPContact")) contact_function=new HPContact();

    Logger(log_thres)<<Name()<<"> Contact function used: "<<cftype<<"\n";

    parse_CF_opts(contopts,contact_function);

    // Now we can set up the contact list

    for (size_t i =0;i<usrcmd.size();++i) {
        parts.clear();
        split(usrcmd[i],parts);
        if (parts[0]==string("file") && parts.size()>=2) {
            ReadMap(parts[1]);
        }

        if (parts[0]=="all") track_all_contacts();
    }

    init();

    return contact_function->init(p);
}

void ContactMap::parse_CF_opts(vector<string> &ops, ContactFunction *f)
{
    vector<string> parts;

    for (size_t iopt=0;iopt<ops.size();++iopt) {
        parts.clear();
        split(ops[iopt],parts);

        if ((parts[0]=="cutoff" or parts[0]=="cut_off") &&parts.size()>=2) {
            f->set_cutoff(strtod(parts[1].c_str(),NULL));
            Logger(log_thres)<<Name()<<"> Cut-off set to "
            <<strtod(parts[1].c_str(),NULL)<<"\n";
        }
    }
}

void ContactMap::avg_fill(int itmp)
{
    Observable::avg_fill(itmp);
    if (gathstat) {
        for (int i=0;i<n_monitored_contacts;++i) {
            if (mainmap[i]) prof[itmp][i]+=1;
        }
        dndt[itmp]++;
    }
}

double ContactMap::evaluate()
{
    mainmap.reset();

    for (int i=0;i<n_monitored_contacts;++i) {
        if ((*contact_function)(contact[i])) mainmap.set(i);
    }

    return (double) mainmap.count();
}

void ContactMap::write_rtkey(Output &op)
{
    Observable::write_rtkey(op);

    if (repSize() ==1) {
        op<<Name()<<": binary state map\n";
    } else {
        for (int i=0;i<repSize();++i)
            op<<Name()<<": binary state map: part "<<i<<"\n";
    }
}

void ContactMap::write_snapshot(Output &op)
{
    op<<"  "<<obsval;

    if (repSize() ==1) {
        op<<"  "<<State();
    } else {
        DivideToGroups();

        for (int i=0;i<repSize();++i)
            op<<"  "<<StateGroup(i);
    }
}

void ContactMap::DivideToGroups()
{
    for (int i=0;i<representation_size;++i) {
        for (int j=0;j<PRFTARGETSIZE;++j) {
            cmap[i][j]=mainmap[i*PRFTARGETSIZE+j];
        }
    }
}

void ContactMap::CombineGroups()
{
    mainmap.reset();

    for (int i=0;i<representation_size;++i) {
        for (int j=0;j<PRFTARGETSIZE;++j) {
            mainmap[i*PRFTARGETSIZE+j]=cmap[i][j];
        }
    }
}

void ContactMap::repSize(int i)
{
    representation_size=std::min(i,PRFMAXCONTACTS/PRFTARGETSIZE);
    cmap.resize(representation_size);
}

void ContactMap::ReadMap(string mapfile)
{
    Logger blog;
    int dummy,i1,i2,lineno=0;
    string line;

    if (!contact.empty()) {
        blog(log_thres)<<Name()<<"> Clearing non-empty contact list!\n";
    }
    contact.clear();
    n_monitored_contacts=0;

    ifstream fin(mapfile.c_str());
    while (getline(fin,line)) {
        istringstream isin(line.c_str());
        ++lineno;
        char first_char;
        isin>>first_char;

        if (line.empty()||first_char=='#') continue;
        else isin.putback(first_char);

        isin>>dummy;isin>>i1;isin>>i2;

        if (isin.fail()) {
            prf::cerr<<"ignoring line "<<lineno
            <<": \""<<isin.str()<<"\"\n";continue;
        }

        ++n_monitored_contacts;
        if (n_monitored_contacts>PRFMAXCONTACTS) {
            prf::cerr<<"Number of contacts to be monitored in file "
            <<mapfile<<" has reached the maximum that can be "
            <<" handled in this build, which is "<<PRFMAXCONTACTS<<"\n"
            <<"If you really need more contacts to be tracked, "
            <<"you should change the value of PRFMAXCONTACTS in "
            <<"model/Observables/ContactMap.hh and recompile PROFASI.\n"
            <<"This run will exit to draw attention to this.\n";
            exit(1);
        }
        contact.push_back(Contact(i1,i2));
    }
    fin.close();

    representation_size=(n_monitored_contacts-1)/PRFTARGETSIZE + 1;
    cmap.resize(representation_size);

    for (int i=0;i<representation_size;++i) cmap[i].reset();

    blog(log_thres)<<Name()<<"> Monitoring "<<n_monitored_contacts
    <<" contacts, using representation with "<<representation_size
    <<" integers.\n";

    blog<<Name()<<"> Start contact list\n";

    for (int i=0;i<n_monitored_contacts;++i) {
        blog<<i<<"  "<<contact[i].first<<"  "<<contact[i].second<<"\n";
    }

    blog<<Name()<<"> End contact list\n";
}

void ContactMap::track_all_contacts()
{
    Logger blog;

    n_monitored_contacts=p->NumberOfLigands();
    if (contact_function->is_symmetric()) {
        n_monitored_contacts*=(n_monitored_contacts-1);
        n_monitored_contacts/=2;
    } else {
        n_monitored_contacts*=n_monitored_contacts;
    }

    if (n_monitored_contacts>=PRFMAXCONTACTS) {
        prf::cerr<<Name()<<"> Tracking of all possible contacts was requested. "
        <<"That number "<<n_monitored_contacts<<" is larger than the maximum "
        <<"that can be handled in this build, which is "<<PRFMAXCONTACTS<<"\n"
        <<"To use more contacts, you should change the value of PRFMAXCONTACTS"
        <<"in model/Observables/ContactMap.hh, and recompile PROFASI with "
        <<"\"make clean\", followed by \"make\". \n\nThis run will exit to "
        <<"draw attention to the above.\n";
        exit(1);
    }

    if (!contact.empty()) {
        blog(log_thres)<<Name()<<"> Clearing non-empty contact list!\n";
    }
    contact.resize(n_monitored_contacts);

    int ictc=0;
    for (int i=0;i<p->NumberOfLigands();++i) {
        if (contact_function->is_symmetric()) {
            for (int j=0;j<i;++j) {
                contact[ictc++]=Contact(i,j);
            }
        } else {
            for (int j=0;j<p->NumberOfLigands();++j) {
                contact[ictc++]=Contact(i,j);
            }
        }
    }

    representation_size=(n_monitored_contacts-1)/PRFTARGETSIZE + 1;
    cmap.resize(representation_size);

    for (int i=0;i<representation_size;++i) cmap[i].reset();

    blog(log_thres)<<Name()<<"> Monitoring "<<n_monitored_contacts
    <<" contacts, using representation with "<<representation_size
    <<" integers.\n";

    if (contact_function->is_symmetric()) {
        blog<<Name()<<"> The contact list consists of all pairs (i,j) of "
                <<"residues with j < i. \n";
    } else {
        blog<<Name()<<"> The contact list consists of all pairs of residues.\n";
    }
}

void ContactMap::WriteMaps()
{
    prf::cout<<"main map : ";

    for (size_t j=0;j<mainmap.size();++j) prf::cout<<(mainmap[j]?'1':'0');

    prf::cout<<"\n";

    for (int i=0;i<representation_size;++i) {
        prf::cout<<"map "<<i<<": ";

        for (size_t j=0;j<cmap[i].size();++j)
            prf::cout<<((bool) cmap[i][j]?'1':'0');

        prf::cout<<"\n";
    }
}

void ContactMap::init()
{
    prof.allocate(ntmp,n_monitored_contacts);
    prof*=0.0;
    dndt.resize(ntmp,0);
}

void ContactMap::avg_reset()
{
    Observable::avg_reset();
    prof*=0.0;
    dndt.assign(ntmp,0);
}

void ContactMap::avg_write(Output &op)
{
    Observable::avg_write(op);
    Output fout((oprefx+Name()+".profile").c_str());
    fout<<"#different rows represent different temperatures\n";
    fout<<"#different columns represent different bonds/contacts\n";

    for (int itmp=0;itmp<ntmp;++itmp) {
        for (int ibnd=0;ibnd<n_monitored_contacts;++ibnd) {
            if (dndt[itmp]) fout<<prof[itmp][ibnd]/dndt[itmp]<<"  ";
            else fout<<"0  ";
        }

        fout<<"\n";
    }

    fout.close();
}
