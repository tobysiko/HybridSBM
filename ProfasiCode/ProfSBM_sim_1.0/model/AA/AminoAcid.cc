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

#include "AminoAcid.hh"

using namespace UnivConstants;
using std::vector;
using std::deque;
using std::pair;
using std::string;

using namespace prf;

bool AminoAcid::initialized=false;
const double AminoAcid::theta[3]={121.7*pi/180,111.0*pi/180,116.6*pi/180};
const double AminoAcid::b[3]={1.33,1.46,1.52};
const double AminoAcid::bNH=1.00;
const double AminoAcid::bCaHa=1.08;
const double AminoAcid::bCaCb=1.53;
const double AminoAcid::bCO=1.23;
const double AminoAcid::phCb=pi-120.9*pi/180;
const double AminoAcid::thCb=110.0*pi/180;
const double AminoAcid::phHa=pi+118.7*pi/180;
const double AminoAcid::thHa=109.0*pi/180;
const double AminoAcid::DPRphHa=pi-118.7*pi/180;
const double AminoAcid::DPRphCb=pi+120.9*pi/180;
const double AminoAcid::phHa2=pi-118.7*pi/180;
FindCoord AminoAcid::locate_Hn;
FindCoord AminoAcid::locate_Ha;
FindCoord AminoAcid::locate_Ha2;
FindCoord AminoAcid::locate_Cb;
FindCoord AminoAcid::locate_Oc;
FindCoord AminoAcid::locate_O2c;
FindCoord AminoAcid::locate_H1n;
FindCoord AminoAcid::locate_H2n;
FindCoord AminoAcid::locate_H3n;
FindCoord AminoAcid::locate_dHa;
FindCoord AminoAcid::locate_dCb;

void AminoAcid::initCommon()
{
    Groups::initGroups();

    if (!initialized) {
        locate_Hn.Initialize(b[0],b[1],bNH,pi-theta[0],
                             (pi-0.5*theta[0]),pi);
        locate_H1n.Initialize(b[0],b[1],bNH,pi-theta[0], (acos(-1.0/3)),pi);
        locate_H2n.Initialize(b[0],b[1],bNH,pi-theta[0],
                              (acos(-1.0/3)),-pi/3);
        locate_H3n.Initialize(b[0],b[1],bNH,pi-theta[0],
                              (acos(-1.0/3)),pi/3);
        locate_Ha.Initialize(b[2],b[1],bCaHa,pi-theta[1],pi-thHa,phHa);
        locate_dHa.Initialize(b[2],b[1],bCaHa,pi-theta[1],
                              pi-thHa,DPRphHa);
        locate_Ha2.Initialize(b[2],b[1],bCaHa,pi-theta[1],pi-thHa,phHa2);
        locate_Cb.Initialize(b[2],b[1],bCaCb,pi-theta[1],pi-thCb,phCb);
        locate_dCb.Initialize(b[2],b[1],bCaCb,pi-theta[1],pi-thCb,DPRphCb);

        locate_Oc.Initialize(b[2],b[0],bCO,pi-theta[2], (pi-0.5*theta[2]),pi);
        locate_O2c.Initialize(b[2],b[0],bCO,pi-theta[2],theta[2]-theta[0],0);

        initialized=true;
    }
}

void AminoAcid::nodes_reconnect()
{
    if (nnodes==0) return;
    int myca=Calpha().UniqueId();
    int n0ca=node[0]->atom(1).UniqueId();
    for (int i=0;i<nnodes;++i) node[i]->AtomOffset(myca-n0ca);
}

void AminoAcid::Initialize() {}

void AminoAcid::Allocate()
{
    ibbc=mygrp->index(" C  ");
    ibbca=mygrp->index(" CA ");
    icb=ibbca+2;
    ioc=mygrp->index(" O  ");
    na=1;

    if (mytype==P||mytype==DPR) nd=0;
    else nd=1;

    if (ntrml) natms+=2;

    if (ctrml) natms+=1;

    Ligand::Allocate();

    int i=mygrp->num_atoms();

    if (ntrml) {    //two more hydrogens if it is N-terminal
        atm[i++].Species(hydrogen);
        atm[i++].Species(hydrogen);
        nd+=2;
    }

    if (ctrml) {atm[i++].Species(oxygen);na++;nhvatm++;}

    //extra charged oxygen at C-terminal
    if (!ntrml) for (int j=0;j<natms;++j) atm[j].UniqueId(j);
    else {
        int aunid=0;
        atm[0].UniqueId(aunid++);

        if (mytype!=P && mytype!=DPR) {
            atm[1].UniqueId(aunid++);
            atm[natms-1].UniqueId(aunid++);
            atm[natms-2].UniqueId(aunid++);
        } else {
            atm[natms-2].UniqueId(aunid++);
            atm[natms-1].UniqueId(aunid++);
        }

        for (int j=ibbca;j<natms-2;++j) atm[j].UniqueId(aunid++);
    }
}

AminoAcid::AminoAcid(OneLetterCode typ) : Ligand(typ)
{
    ntrml=ctrml=cis=false;
    nsatm=natms- ((mytype==P||mytype==DPR) ?5:6);
    theBB=NULL;
    seqser=0;
    BBloc=chser=-1;
}

AminoAcid::AminoAcid(const AminoAcid &gac) : Ligand(gac), ntrml(gac.ntrml),
        ctrml(gac.ctrml), nsatm(gac.nsatm), ibbca(gac.ibbca), ibbc(gac.ibbc),
        ioc(gac.ioc) , icb(gac.icb), theBB(gac.theBB),
        nnodes(gac.nnodes), BBloc(gac.BBloc), chirality(gac.chirality) {}

AminoAcid & AminoAcid::operator=(const AminoAcid &gac)
{
    if (this!=&gac) {
        Ligand::operator=(gac);
        ntrml=gac.ntrml;
        ctrml=gac.ctrml;
        nsatm=gac.nsatm;
        ibbca=gac.ibbca;
        ibbc=gac.ibbc;
        ioc=gac.ioc;
        icb=gac.icb;
        theBB=gac.theBB;
        BBloc=gac.BBloc;
        nnodes=gac.nnodes;
        chirality=gac.chirality;
    }

    return *this;
}

AminoAcid::~AminoAcid() {}

void AminoAcid::Reconstruct()
{
    //Assumes that the backbone is in place. Make sure back bone reconstruction
    //is called before this routine.
    //Proline does not depend on this function. When this function
    //is reached, the amino acid is definitely not D- or L-Proline.
    if (hasNTerminal()) {
        atm[1].Pos(Nitrogen().Pos()+
                    locate_H1n(theBB->bond(BBloc),theBB->bond(BBloc+1)));
        atm[natms-1].Pos(Nitrogen().Pos() +
                    locate_H2n(theBB->bond(BBloc),theBB->bond(BBloc+1)));
        atm[natms-2].Pos(Nitrogen().Pos() +
                    locate_H3n(theBB->bond(BBloc),theBB->bond(BBloc+1)));
    } else atm[1].Pos(Nitrogen().Pos() +
                    locate_Hn(theBB->bond(BBloc),theBB->bond(BBloc+1)));

    atm[ibbca+1].Pos(Calpha().Pos() +
        locate_Ha(theBB->bond(BBloc+2),theBB->bond(BBloc+1)));

    atm[ibbc+1].Pos(Cprime().Pos() +
        locate_Oc(theBB->bond(BBloc+2),theBB->bond(BBloc+3)));

    if (hasCTerminal()) {
        atm[ibbc+2].Pos(Cprime().Pos() +
            locate_O2c(theBB->bond(BBloc+2),theBB->bond(BBloc+3)));
    }

    if (mytype!=G) {
        atm[ibbca+2].Pos(Calpha().Pos() +
                    locate_Cb(theBB->bond(BBloc+2),theBB->bond(BBloc+1)));
        node[0]->ReCreate();
    } else
        atm[ibbca+2].Pos(Calpha().Pos()+
            locate_Ha2(theBB->bond(BBloc+2),theBB->bond(BBloc+1)));
}

int AminoAcid::UpdateSideChainDof(int rdof) {return rdof;}

Atom * AminoAcid::labeled_atom(string alabel)
{
    Logger blog;

    if (hasNTerminal()) {
        if (alabel==string("2H  ")) return &atm[natms-1];

        if (alabel==string("3H  ")) {
            if (mytype!=P&&mytype!=DPR) return &atm[natms-2];
            else return NULL;
        }

        if (alabel==string("1H  ")) {
            if (mytype!=P&&mytype!=DPR) return &atm[1];
            else return &atm[natms-2];
        }
    }

    if (hasCTerminal() &&alabel==string(" OXT")) return &atm[ioc+1];

    int indx=mygrp->index(alabel);

    if (indx>=0) return &atm[indx];

    blog(500) <<"Warning: could not locate atom for label "<<alabel
    <<" on amino acid "<<Name()
    <<"("<<seqser<<"). NULL pointer returned. \n";

    return NULL;
}

string AminoAcid::label_of(int i)
{
    if (ntrml) {
        if (i==1 && ibbca==2) return string("1H  ");
        else if (i==natms-2) return ibbca==2?"3H  ":"1H  ";
        else if (i==natms-1) return "2H  ";
    } else if (ctrml) {
        if (i==natms-1) return string(" OXT");
    }

    return mygrp->label(i);
}

Atom & AminoAcid::at(string alabel)
{
    if (hasNTerminal()) {
        if (alabel==string("2H  ")) return atm[natms-1];

        if (alabel==string("3H  ")) {
            if (mytype!=P&&mytype!=DPR) return atm[natms-2];
            else return atm[natms-1]; //Should actually be NULL, but ...
        }

        if (alabel==string("1H  ")) {
            if (mytype!=P&&mytype!=DPR) return atm[1];
            else return atm[natms-2];
        }
    }

    if (hasCTerminal() &&alabel==string(" OXT")) return atm[ioc+1];

    return Ligand::at(alabel);
}

string AminoAcid::USCAtomDescr(int iscatm)
{
    switch (atm[iscatm+ibbca+1].Species()) {

        case hydrogen:
            return string(" HX ");
        case carbon:
            return string(" CX ");
        case nitrogen:
            return string(" NX ");
        case oxygen:
            return string(" OX ");
        case sulfur:
            return string(" SX ");
        default:
            return string(" ERR");
    };
}

void AminoAcid::Write()
{
    prf::cout <<Name();
    prf::cout <<" ("<<CharCode() <<") ";

    if (hasNTerminal()) prf::cout <<" (N-terminal) ";

    if (hasCTerminal()) prf::cout <<" (C-terminal) ";

    prf::cout <<"\n";

    for (int i=0;i<NumberOfAtoms();++i) atm[i].Write();
}

void AminoAcid::WritePDB(int &atindx, char ch_id,int aaindx,FILE *fp)
{
    WritePDBline(atindx,ch_id,aaindx,0,fp);

    if (ibbca!=1) WritePDBline(atindx,ch_id,aaindx,1,fp);

    int i=0;

    if (ntrml) {
        for (i=atm.size()-1;i>=mygrp->num_atoms();--i)
            WritePDBline(atindx,ch_id,aaindx,i,fp);
    }

    for (i=ibbca;i<mygrp->num_atoms();++i)
        WritePDBline(atindx,ch_id,aaindx,i,fp);

    if (ctrml) WritePDBline(atindx,ch_id,aaindx,i,fp);
}

void AminoAcid::WritePDB2(int &atindx, char ch_id,int aaindx,FILE *fp)
{
    WritePDBline(atindx,ch_id,aaindx,0,fp);
    WritePDBline(atindx,ch_id,aaindx,ibbca,fp);
    WritePDBline(atindx,ch_id,aaindx,ibbc,fp);
    WritePDBline(atindx,ch_id,aaindx,ioc,fp);

    for (int i=0;i<numSideChainAtoms();++i) {
        if (atm[ibbca+2+i].Species() !=hydrogen) {
            WritePDBline(atindx,ch_id,aaindx,ibbca+2+i,fp);
        }
    }

    if (hasCTerminal()) WritePDBline(atindx,ch_id,aaindx,ioc+1,fp);

    for (int i=0;i<NumberOfAtoms();++i) {
        if (atm[i].Species() ==hydrogen)
            WritePDBline(atindx,ch_id,aaindx,i,fp);
    }
}

int AminoAcid::rotDof_assign_and_reconstruct(int i, double mgd, int &a0, int &a1)
{
    node[i]->AssignPhi(mgd);
    node[i]->MobileAtoms(a0,a1);
    return 1;
}

int AminoAcid::rotDof_assign(int il, double mgd)
{
    node[il]->RevertPhi(mgd);
    return 1;
}

double AminoAcid::get_rotDof(int il)
{
    return node[il]->Phi();
}

int AminoAcid::ROTDOF(int il, double mgd, int &a0, int &a1)
{

    int j;
    node[j= (il-rtOffset)]->AssignPhi(mgd);
    node[j]->MobileAtoms(a0,a1);
    return 1;
}

int AminoAcid::ROTDOFr(int il, double mgd)
{
    node[(il-rtOffset)]->RevertPhi(mgd);
    return 1;
}

double AminoAcid::ADOF(int il)
{
    return node[il-rtOffset]->Phi();
}

void AminoAcid::LocPairsatRTdof(int i, deque<pair<int,int> > & lcp)
{
    node[(i-rtOffset)]->LocPairs(lcp);
}

double AminoAcid::get_coord(int i)
{
    switch (i) {
        case 0: return Phi();
        case 1: return Psi();
        case 2: return Omega();
        default: return get_rotDof(i-3);
    };
}

void AminoAcid::set_coord(int i, double x)
{
    switch (i) {
        case 0: Phi(x); break;
        case 1: Psi(x); break;
        case 2: Omega(x); break;
        default: rotDof_assign(i-3,x);
    };
}

int AminoAcid::n_coord() const {
    return 3+n_rotDof();
}

int AminoAcid::n_dof() const {
    return (mytype==P or mytype==DPR)?1:2+n_rotDof();
}

double AminoAcid::get_dof(int i)
{
    if (mytype==P or mytype==DPR) return Psi();
    switch (i) {
        case 0: return Phi();
        case 1: return Psi();
        default: return get_rotDof(i-2);
    };
}

void AminoAcid::set_dof(int i, double x)
{
    if (mytype==P or mytype==DPR) Psi(x);
    else {
        switch (i) {
        case 0: Phi(x); break;
        case 1: Psi(x); break;
        default: rotDof_assign(i-2,x);
        };
    }
}

void AminoAcid::WriteNodes()
{
    int ist,ind;

    for (unsigned int i=0;i<node.size();++i) {
        if (node[i]!=NULL) {
            node[i]->MobileAtoms(ist,ind);
            prf::cout <<"node "<<i<<"("<<ist<<", "<<ind<<")\n";
        } else {
            prf::cout <<"NULL Node\n";
        }
    }
}

void AminoAcid::WriteNodes(string &nds)
{
    int ist,ind;
    char tmp[80];

    for (unsigned int i=0;i<node.size();++i) {
        if (node[i]!=NULL) {
            node[i]->MobileAtoms(ist,ind);
            sprintf(tmp,"node %d %s (%d, %d)\n",i,node[i]->Name().c_str(),ist,ind);
        } else {
            sprintf(tmp,"node %d is NULL\n",i);
        }

        nds+=tmp;
    }
}

void AminoAcid::Donor(int i, Dipole &dp)
{
    if (nd>1) {
        switch (i%nd) {

            case 0:
                dp.SetAtoms(at(" N  ").UniqueId(),at("1H  ").UniqueId());
                break;

            case 1:
                dp.SetAtoms(at(" N  ").UniqueId(),at("2H  ").UniqueId());
                break;

            case 2:
                dp.SetAtoms(at(" N  ").UniqueId(),at("3H  ").UniqueId());
                break;

            default:
                break;
        }
    } else if (nd==1)
        dp.SetAtoms(at(" N  ").UniqueId(),at(" H  ").UniqueId());
}

void AminoAcid::Acceptor(int i, Dipole &dp)
{
    if (na>0) {
        switch (i%na) {

            case 0:
                dp.SetAtoms(atm[ibbc+1].UniqueId(),
                            atm[ibbc].UniqueId());
                break;

            case 1:
                dp.SetAtoms(atm[ibbc+2].UniqueId(),
                            atm[ibbc].UniqueId());
                break;

            default:
                break;
        }
    }
}

void AminoAcid::ExportConnections(ConnectionsMatrix &aa)
{
    if (mytype==P||mytype==DPR) {
        for (int i=0;i<13;++i) {
            for (int j=i;j<13;++j) {
                aa.set_connection(atm[0].UniqueId() +i,atm[0].UniqueId() +j);
            }
        }
    } else if (mytype!=G) {
        node[0]->ExportConnections(aa);
    }

    if (mytype==Y) {
        //Temporary solution until I find a better way to export these
        //connections. These atoms are really far in terms of bonds,
        //and yet have a fixed distance in the model.
        aa.set_connection(atm[ibbca+2].UniqueId(),atm[ibbca+16].UniqueId());
        aa.set_connection(atm[ibbca+5].UniqueId(),atm[ibbca+16].UniqueId());
    }
}

void AminoAcid::BuildConnections()
{
    if (mytype!=G) {
        if (node[0]!=NULL) {
            node[0]->set_bases(Nitrogen(), Hca(),Cprime());
            node[0]->BuildConnections();
        }
    }
}

void AminoAcid::calcPhiPsi(double &phv, double &psv) const
{
    phv=theBB->bond(BBloc+1).torsion(theBB->bond(BBloc),theBB->bond(BBloc+2));
    psv=theBB->bond(BBloc+2).torsion(theBB->bond(BBloc+1),
                                    theBB->bond(BBloc+3));
}

bool AminoAcid::is_helical() const
{
    double phv=Phi(),psv=Psi();
    return ((phv>helix_phimin)&&(phv<helix_phimax)
            &&(psv>helix_psimin)&&(psv<helix_psimax));
}

bool AminoAcid::is_strand() const
{
    double phv=Phi(),psv=Psi();
    return ((phv>sheet_phimin)&&(phv<sheet_phimax)
            &&(psv>sheet_psimin)&&(psv<sheet_psimax));
}
