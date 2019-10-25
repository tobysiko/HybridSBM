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

#include "../Aux/Constants.hh"
#include "Protein.hh"
#include "../Aux/SeqBuilder.hh"
#include <fstream>
#include <sstream>

using namespace prf;

using namespace UnivConstants;
using std::ifstream;
using std::string;
using std::vector;
using std::deque;
using std::pair;
using std::make_pair;

Protein::Protein() : bb(0)
{
    NtLigand=CtLigand=NULL;
    clear();
}

Protein::Protein(const Protein & pptd) : bb(0)
{
    NtLigand=CtLigand=NULL;
    charged_NT=pptd.charged_NT;
    charged_CT=pptd.charged_CT;
    set_sequence(pptd.seq);
    unid=pptd.unid;
    setAtomOffset(pptd.atoffset);
    consi=pptd.consi;
    std::vector<double> tmpdof;
    pptd.get_dof(tmpdof);
    set_dof(tmpdof);
}

Protein & Protein::operator= (const Protein & pptd)
{
    if (this!=&pptd) {
        if (Sequence()!=pptd.Sequence()) {
            clear();
            charged_NT=pptd.charged_NT;
            charged_CT=pptd.charged_CT;
            set_sequence(pptd.seq);
        }
        unid=pptd.unid;
        setAtomOffset(pptd.atoffset);
        consi=pptd.consi;
        std::vector<double> tmpdof;
        pptd.get_dof(tmpdof);
        set_dof(tmpdof);
    }

    return *this;
}

Protein::Protein(string aaseq) :bb(0)
{
    NtLigand=CtLigand=NULL;
    charged_NT=charged_CT=true;
    Allocate(0,0,aaseq);
}

Protein::~Protein()
{
    clear();
}

void Protein::clear()
{
    if (NtLigand!=NULL) delete NtLigand;

    if (!AAv.empty()) {
        for (unsigned int i=0;i<AAv.size();++i)
            if (AAv[i]!=NULL) delete AAv[i];
    }

    if (CtLigand!=NULL) delete CtLigand;
    AAv.clear();
    lg.clear();
    bb.clear();
    seq.clear();
    rtloc.clear();
    cisres.clear();
    NtLigand=CtLigand=NULL;
    unid=0;
    atbg=atnd=lgbg=lgnd=naa=nrtval=natms=atoffset=0;
    charged_NT=charged_CT=true;
    consi=false;
}

std::string Protein::Sequence() const
{
    SeqBuild sbuild;
    return sbuild.make_string(seq,' ',' ');
}

void Protein::setUniqueId(int i)
{
    unid=i;
    for (size_t j=0;j<lg.size();++j) lg[j]->LocatedOn(unid);
}

void Protein::setAtomOffset(int i)
{
    atoffset=i;
    for (size_t j=0;j<lg.size();++j) {
        lg[j]->atomOffset(i);
        i+=lg[j]->NumberOfAtoms();
    }
    atbg=atoffset;
    atnd=atbg+natms;
    for (size_t j=0;j<AAv.size();++j) AAv[j]->nodes_reconnect();
    if (NtLigand!=NULL) {
        NtLigand->nodes_reconnect(AAv[0]->Calpha(),AAv[0]->Nitrogen());
    }

    if (CtLigand!=NULL) {
        CtLigand->nodes_reconnect(AAv[naa-1]->Calpha(),AAv[naa-1]->Cprime());
    }

    bb.set_atom_offset(AAv[0]->Nitrogen().UniqueId());
}

int Protein::set_sequence(const std::vector<prf::OneLetterCode> &gseq)
{
    if (gseq.empty()) return -1;

    if (Groups::grp[gseq[0]].isEG()) add_NTLigand(gseq[0]);
    else if (!charged_NT) add_NTLigand(VOIDEG);
    for (size_t i=0;i<gseq.size();++i) {
        if (Groups::grp[gseq[i]].isAA()) add_aminoacid(gseq[i]);
    }
    if (Groups::grp[gseq.back()].isEG()) add_CTLigand(gseq.back());
    else if (!charged_CT) add_CTLigand(VOIDEG);

    connect_residues();
    init_atoms();
    create_backbone();
    init_dof();
    for (size_t i=0;i<lg.size();++i) lg[i]->Initialize();
    BuildConnections();
    return 1;
}

void Protein::connect_residues()
{
    lg.clear();
    if (NtLigand!=NULL) {
        lg.push_back(NtLigand);
    }
    for (size_t i=0;i<AAv.size();++i) {
        lg.push_back(AAv[i]);
        AAv[i]->hasNTerminal(false);
        AAv[i]->hasCTerminal(false);
    }
    if (CtLigand!=NULL) {
        lg.push_back(CtLigand);
    }

    if (NtLigand==NULL && !AAv.empty()) AAv.front()->hasNTerminal(true);
    if (CtLigand==NULL && !AAv.empty()) AAv.back()->hasCTerminal(true);

    for (size_t i=0;i<lg.size();++i) {
        if (i>0) lg[i]->LeftConnection(lg[i-1]);
        if (i<lg.size()-1) lg[i]->RightConnection(lg[i+1]);
        lg[i]->SeqSerial(i);
    }
}

void Protein::create_backbone()
{
    bb.clear();
    bb.numAminoAcids(AAv.size());
    for (int i=0;i<naa;++i) bb.registerAA(i,AAv[i]);
    if (NtLigand!=NULL) {
        NtLigand->Connect_to_AA(AAv[0]);
        NtLigand->SetRefv(bb.Bond(1),bb.Bond(0));
        NtLigand->SetRefa(AAv[0]->Calpha(),AAv[0]->Nitrogen());
    }

    if (CtLigand!=NULL) {
        CtLigand->Connect_to_AA(AAv[naa-1]);
        CtLigand->SetRefv(bb.Bond(3*naa-1),bb.Bond(3*naa));
        CtLigand->SetRefa(AAv[naa-1]->Calpha(),AAv[naa-1]->Cprime());
    }
}

void Protein::init_dof()
{
    Logger blog(10);
    nrtval=0;
    for (size_t i=0;i<lg.size();++i) {
        lg[i]->rotOffset(nrtval);
        nrtval+=lg[i]->n_rotDof();
    }

    rtloc.resize(nrtval,0);

    int j=0;

    for (size_t i=0;i< lg.size();++i) {
        int k=lg[i]->n_rotDof();
        while (k-->0) rtloc[j++]=i;
    }

    if (nrtval!=j) {
        prf::cerr<<"mismatch in number of rotamer degrees of freedom counted "
        <<"and those that were allocated. \n";

    }
    for (size_t i=0;i<cisres.size();++i) {
        if (cisres[i]>=0 && cisres[i]<(naa-1)) { //last residue cannot be cis
            blog(5)<<"Protein> setting residue "<<cisres[i]<<" to cis\n";
            AAv[cisres[i]]->setCis();
        } else {
            prf::cerr<<"Protein> ignoring illegal cis-request for residue"
            <<cisres[i]<<" in chain "<<unid<<"\n";
        }
    }
}

void Protein::init_atoms()
{
    natms=0;
    for (size_t i=0;i<lg.size();++i) {
        lg[i]->Allocate();
        lg[i]->atomOffset(natms);
        natms+=lg[i]->NumberOfAtoms();
    }
    atbg=atoffset;
    atnd=atbg+natms;
}

int Protein::add_aminoacid(prf::OneLetterCode cd)
{
    if (not Groups::grp[cd].isAA()) return 0;
    AminoAcid *ac=new_AA_object(cd);
    ac->LocatedOn(unid);
    AAv.push_back(ac);
    ++naa;
    seq.push_back(cd);
    return 1;
}

int Protein::add_NTLigand(prf::OneLetterCode cd)
{
    if (not Groups::grp[cd].isEG()) return 0;
    if (NtLigand!=NULL) delete NtLigand;
    NtLigand=new_EG_object(cd);
    NtLigand->LocatedOn(unid);
    if (seq.empty() or !Groups::grp[seq[0]].isEG()) {
        std::vector<prf::OneLetterCode> tmps(seq.size()+1,cd);
        std::copy(seq.begin(),seq.end(),tmps.end());
        seq=tmps;
    } else seq[0]=cd;
    return 1;
}

int Protein::add_CTLigand(prf::OneLetterCode cd)
{
    if (not Groups::grp[cd].isEG()) return 0;
    if (CtLigand!=NULL) delete CtLigand;
    CtLigand=new_EG_object(cd);
    CtLigand->LocatedOn(unid);
    if (seq.empty() or !Groups::grp[seq.back()].isEG()) seq.push_back(cd);
    else seq[seq.size()-1]=cd;
    return 1;
}

void Protein::Allocate(int ppid,int atom_offset,string sq)
{
    vector<OneLetterCode> codes;
    SeqBuild sbuild;
    sbuild.parse(sq,codes);
    Allocate(ppid,atom_offset,codes);
}

void Protein::Allocate(int ppid,int atom_offset,vector<OneLetterCode> fulseq)
{
    if (fulseq.empty()) return;
    bool b1=charged_NT,b2=charged_CT;
    std::vector<int> cisresb=cisres;
    clear();
    charged_NT=b1;charged_CT=b2;
    cisres=cisresb;
    if (not good_sequence(fulseq)) return;
    unid=ppid;
    set_sequence(fulseq);
    setAtomOffset(atom_offset);
    consi=true;
}

bool Protein::good_sequence(vector<OneLetterCode> sq)
{
    bool ans=true;
    for (size_t i=0;i<sq.size();++i) {
        if (i==0 or i==(sq.size()-1)) {
            ans=ans&&(Groups::grp[sq[i]].isEG() or Groups::grp[sq[i]].isAA());
        } else {
            ans=ans&&Groups::grp[sq[i]].isAA();
        }
    }
    return ans;
}
void Protein::reconstruct()
{
    bb.forwardReconstruct(0,3*naa+1);
    //bb.reverseReconstruct(3*naa+1,0);

    for (size_t i=0;i<lg.size();++i) lg[i]->Reconstruct();
}

void Protein::reconstruct(int idir, int iaa)
{
    if (idir==0) {
        bb.forwardReconstruct(3*iaa,3*naa+1);

        for (int i=iaa;i<naa;++i) AAv[i]->Reconstruct();

        if (CtLigand!=NULL) CtLigand->Reconstruct();

        if (iaa==0 && NtLigand!=NULL) NtLigand->Reconstruct();
    } else {
        bb.reverseReconstruct(3*iaa+3,0);

        for (int i=iaa;i>=0;--i) AAv[i]->Reconstruct();

        if (NtLigand!=NULL) NtLigand->Reconstruct();

        if (iaa== (naa-1) && CtLigand!=NULL) CtLigand->Reconstruct();
    }
}

void Protein::trivial_init()
{
    for (int i=0;i<numBBdof();++i) BBdof(i,0.0);

    for (int i=0;i<numRTdof();++i) RTdof(i,0.0);

    reconstruct();
}

void Protein::stretched_init()
{
    for (int i=0;i<numBBdof();++i) {
        if (backbone()->DOFid(i) %3==0) {
            BBdof(i,-pi);
        } else BBdof(i,pi);
    }

    for (int i=0;i<numRTdof();++i) RTdof(i,pi);

    reconstruct();
}

void Protein::helical_init()
{
    backbone()->SetToAlphaHelix();

    for (int i=0;i<numRTdof();++i) RTdof(i,pi);

    reconstruct();
}

void Protein::randomize(RandomNumberBase *rangen)
{
    for (int i=0;i<numBBdof();++i) BBdof(i,twoPi*(rangen->shoot()));

    for (int i=0;i<numRTdof();++i) RTdof(i,twoPi*(rangen->shoot()));

    reconstruct();
}

void Protein::Write()
{
    prf::cout<<"Writing peptide "<< Sequence() <<"\n";

    if (NtLigand!=NULL) NtLigand->Write();

    prf::cout<<"Amino Acid table \n";

    for (int i=0;i<naa;++i) AAv[i]->Write();

    if (CtLigand!=NULL) CtLigand->Write();
}

void Protein::WritePDB(int &atindx,FILE *fp,char ch_id, int &rsindx)
{
    int aaindx=rsindx;
    size_t j=0;

    for (j=0;j<lg.size();++j) {
        if (lg[j]) lg[j]->WritePDB(atindx,ch_id,++aaindx,fp);
    }

    fprintf(fp,"TER    %4u      %s %c%4u\n",
            ++atindx,lg.back()->TLC().c_str(),ch_id,aaindx);
    rsindx=aaindx;
}

void Protein::WritePDB2(int &atindx,FILE *fp,char ch_id, int &rsindx)
{
    int aaindx=rsindx;
    size_t j=0;

    for (j=0;j<lg.size();++j) {
        if (lg[j]) {
            lg[j]->WritePDB2(atindx,ch_id,++aaindx,fp);
        }
    }

    fprintf(fp,"TER    %4u      %s %c%4u\n",
            ++atindx,lg.back()->TLC().c_str(),ch_id,aaindx);
    rsindx=aaindx;
}

void Protein::ExportConnections(ConnectionsMatrix &aa)
{
    int i,j,n;
    Atom *at1=NULL,*at2=NULL;

    for (i=0;i<natms;i++) aa.set_connection(i,i);

    //All intra-residue connections of up to 3 covalent bonds.
    string lst[]={" N  "," CA "," HA "," CB ","1HA ","2HA "," C  ",
                  " O  "," OXT"," H  ","1H  ","2H  ","3H  "
                 };

    //x from list vs N from next residue
    string lstxN[]={" N  "," HA ","1HA ","2HA "," CB "," CA "};

    //CA vs x from list in next residue
    string lstCax[]={" CA "," H  "};

    //O vs x from list in next residue
    string lstOx[]={" N  "," H  "," CA "};

    //C vs x from list in next residue
    string lstCx[]={" N  "," H  "," CA "," HA ","1HA ",
                    "2HA "," CB "," C  "
                   };

    for (n=0;n<naa;++n) {
        //Connections within residue
        for (i=0;i<7;++i) {
            for (j=0;j<13;++j)
                if ((at1=AA(n)->labeled_atom(lst[i])) !=NULL &&
                    (at2=AA(n)->labeled_atom(lst[j])) !=NULL)
                    aa.set_connection(at1->UniqueId(),at2->UniqueId());
        }

        for (i=10;i<13;++i) {
            for (j=10;j<13;++j)
                if ((at1=AA(n)->labeled_atom(lst[i])) !=NULL &&
                    (at2=AA(n)->labeled_atom(lst[j])) !=NULL)
                    aa.set_connection(at1->UniqueId(),at2->UniqueId());
        }

        at1=AA(n)->labeled_atom(" OXT");

        if (at1) aa.set_connection(AA(n)->labeled_atom(" O  ")->UniqueId(),
                                       at1->UniqueId());

        //connections to the next residue
        if (n<naa-1) {
            j=AA(n+1)->labeled_atom(" N  ")->UniqueId();

            for (i=0;i<6;++i)
                if ((at1=AA(n)->labeled_atom(lstxN[i])) !=NULL)
                    aa.set_connection(at1->UniqueId(),j);

            i=AA(n)->labeled_atom(" CA ")->UniqueId();

            for (j=0;j<2;++j)
                if ((at2=AA(n+1)->labeled_atom(lstCax[j])) !=NULL)
                    aa.set_connection(i,at2->UniqueId());

            i=AA(n)->labeled_atom(" O  ")->UniqueId();

            for (j=0;j<3;++j)
                if ((at2=AA(n+1)->labeled_atom(lstOx[j])) !=NULL)
                    aa.set_connection(i,at2->UniqueId());

            i=AA(n)->labeled_atom(" C  ")->UniqueId();

            for (j=0;j<8;++j)
                if ((at2=AA(n+1)->labeled_atom(lstCx[j])) !=NULL)
                    aa.set_connection(i,at2->UniqueId());
        }
    }

    /* distances fixed because of fixed Pro phi angle */
    for (i=1;i<naa;i++) {
        if (seq[i]==P||seq[i]==DPR) {
            for (j=0;j<13;++j) {
                aa.set_connection(AA(i-1)->Calpha().UniqueId(),iN(i) +j);
                aa.set_connection(AA(i-1)->Cprime().UniqueId(),iN(i) +j);
                aa.set_connection(AA(i-1)->Oc().UniqueId(),iN(i) +j);
            }
        }
    }

    ExportConnectionssc(aa);
}

void Protein::ExportConnectionssc(ConnectionsMatrix &aa)
{
    for (size_t i=0;i<lg.size();++i) lg[i]->ExportConnections(aa);
}

void Protein::EnforceBC()
{
    bb.FixAngles();
    AtomCoordinates::EnforceBC(atbg,atnd);
}

void Protein::WriteConf(FILE *fp)
{
    double xx=0,yy=0,zz=0;

    for (int i=0;i<3;++i) {
        xx=bb.memberAtom(i).Pos().x();
        yy=bb.memberAtom(i).Pos().y();
        zz=bb.memberAtom(i).Pos().z();
        fwrite(&xx,sizeof(double),1,fp);
        fwrite(&yy,sizeof(double),1,fp);
        fwrite(&zz,sizeof(double),1,fp);
    }

    for (int i=0;i<bb.numDOF();++i) {
        xx=bb.DOF(i);
        fwrite(&xx,sizeof(double),1,fp);
    }

    for (int i=0;i<nrtval;++i) {
        xx=RotDOF(i);
        fwrite(&xx,sizeof(double),1,fp);
    }
}

void Protein::ReadConf(FILE *fp)
{
    double xx=0,yy=0,zz=0;
    Vector3 vv;

    for (int i=0;i<3;++i) {
        fread(&xx,sizeof(double),1,fp);
        fread(&yy,sizeof(double),1,fp);
        fread(&zz,sizeof(double),1,fp);
        vv=Vector3(xx,yy,zz);
        bb.memberAtom(i).Pos(vv);
    }

    for (int i=0;i<bb.numDOF();++i) {
        fread(&xx,sizeof(double),1,fp);
        bb.DOF(i,xx);
    }

    for (int i=0;i<nrtval;++i) {
        fread(&xx,sizeof(double),1,fp);
        RotDOFr(i,xx);
    }
}

void Protein::WriteConf_text(FILE *fp)
{
    double xx=0,yy=0,zz=0;

    for (int i=0;i<3;++i) {
        xx=bb.memberAtom(i).Pos().x();
        yy=bb.memberAtom(i).Pos().y();
        zz=bb.memberAtom(i).Pos().z();
        fprintf(fp,"%.16f\n",xx);
        fprintf(fp,"%.16f\n",yy);
        fprintf(fp,"%.16f\n",zz);
    }

    for (int i=0;i<bb.numDOF();++i) {
        xx=bb.DOF(i);
        fprintf(fp,"%.16f\n",xx);
    }

    for (int i=0;i<nrtval;++i) {
        xx=RotDOF(i);
        fprintf(fp,"%.16f\n",xx);
    }
}

void Protein::ReadConf_text(FILE *fp)
{
    double xx=0,yy=0,zz=0;
    Vector3 vv;

    for (int i=0;i<3;++i) {
        fscanf(fp,"%lf",&xx);
        fscanf(fp,"%lf",&yy);
        fscanf(fp,"%lf",&zz);
        vv=Vector3(xx,yy,zz);
        bb.memberAtom(i).Pos(vv);
    }

    for (int i=0;i<bb.numDOF();++i) {
        fscanf(fp,"%lf",&xx);
        bb.DOF(i,xx);
    }

    for (int i=0;i<nrtval;++i) {
        fscanf(fp,"%lf",&xx);
        RotDOFr(i,xx);
    }
}

string Protein::ConfSignature()
{
    std::ostringstream sout;
    sout<<"sequence "<<Sequence()<<"\n";
    sout<<"double "<<(9+numBBdof()+numRTdof())
    <<" "<<sizeof(double)/sizeof(char)<<" "
    <<(9+numBBdof()+numRTdof())*sizeof(double)/sizeof(char)<<"\n";
    return sout.str();
}

void Protein::LocPairsBBdof(int i, deque<pair<int,int> > & lcp)
{
    int iaa, iphi;
    iaa=bb.DOFloc(i);
    iphi=bb.DOFid(i);
    lcp.clear();

    if ((AAv[iaa]->OLC() ==P||
         AAv[iaa]->OLC() ==DPR) && iphi%3==0) return;

    if (iphi%3==2) return;

    if (iphi%3==0) {
        string hlbls[]={" H  ","1H  ","2H  ","3H  "};

        for (int hi=AAv[iaa]->hasNTerminal() ?1:0;hi<4;++hi) {
            Atom *aat=AAv[iaa]->labeled_atom(hlbls[hi]);

            if (aat==NULL) continue;

            int j=aat->UniqueId();

            lcp.push_back(make_pair(j,AAv[iaa]->Hca().UniqueId()));

            lcp.push_back(make_pair(j,AAv[iaa]->Cbeta().UniqueId()));

            lcp.push_back(make_pair(j,AAv[iaa]->Cprime().UniqueId()));
        }

        if (iaa==0) {
            if (NtLigand!=NULL and NtLigand->OLC()!=VOIDEG) {
                lcp.push_back(make_pair(NtLigand->ATOM(0).UniqueId(),
                                        AAv[iaa]->Hca().UniqueId()));
                lcp.push_back(make_pair(NtLigand->ATOM(0).UniqueId(),
                                        AAv[iaa]->Cbeta().UniqueId()));
                lcp.push_back(make_pair(NtLigand->ATOM(0).UniqueId(),
                                        AAv[iaa]->Cprime().UniqueId()));
            }
        } else {
            int j=AAv[iaa-1]->Cprime().UniqueId();
            lcp.push_back(make_pair(j,AAv[iaa]->Hca().UniqueId()));
            lcp.push_back(make_pair(j,AAv[iaa]->Cbeta().UniqueId()));
            lcp.push_back(make_pair(j,AAv[iaa]->Cprime().UniqueId()));
        }
    } else {
        int j=AAv[iaa]->Nitrogen().UniqueId();
        lcp.push_back(make_pair(j,AAv[iaa]->Oc().UniqueId()));
        if ((iaa<naa-1) or AAv[iaa]->hasCTerminal() or CtLigand->OLC()!=VOIDEG)
            lcp.push_back(make_pair(j,AAv[iaa]->Oc().UniqueId() +1));

        j=AAv[iaa]->Hca().UniqueId();
        lcp.push_back(make_pair(j,AAv[iaa]->Oc().UniqueId()));
        if ((iaa<naa-1) or AAv[iaa]->hasCTerminal() or CtLigand->OLC()!=VOIDEG)
            lcp.push_back(make_pair(j,AAv[iaa]->Oc().UniqueId() +1));

        j=AAv[iaa]->Cbeta().UniqueId();
        lcp.push_back(make_pair(j,AAv[iaa]->Oc().UniqueId()));
        if ((iaa<naa-1) or AAv[iaa]->hasCTerminal() or CtLigand->OLC()!=VOIDEG)
            lcp.push_back(make_pair(j,AAv[iaa]->Oc().UniqueId() +1));
    }
}

void Protein::LocPairsRTdof(int i, deque<pair<int,int> > & lcp)
{
    Logger blog;
    int ilg;
    ilg=rtloc[i];
    blog(100) <<"Protein: rt dof "<<i<<" is in ligand "<<ilg<<"\n";
    lcp.clear();

    if (lg[ilg]->OLC() ==P||lg[ilg]->OLC() ==DPR) return;

    lg[ilg]->LocPairsatRTdof(i,lcp);
}

void Protein::BuildConnections()
{
    for (size_t i=0;i<lg.size();++i) {
        lg[i]->BuildConnections();
    }
}

double Protein::get_dof(int i)
{
    if (i<9) {
        int j=i/3,k=i%3;
        return backbone_atom(j).Position()[k];
    }
    i-=9;
    if (i<numBBdof()) return BBdof(i);
    i-=numBBdof();
    if (i<numRTdof()) return RTdof(i);
    return 0;
}

void Protein::set_dof(int i, double vl)
{
    if (i<9) {
        int j=i/3,k=i%3;
        switch (k) {
          case 0: backbone_atom(j).Pos().x(vl); break;
          case 1: backbone_atom(j).Pos().y(vl); break;
          case 2: backbone_atom(j).Pos().z(vl); break;
          default: ;
        };
    } else if (i<(9+numBBdof())) BBdof(i-9,vl);
    else if (i<(9+numBBdof()+numRTdof())) RTdof(i-9-numBBdof(),vl);
}

double Protein::get_coord(int i)
{
    if (i<9) {
        int j=i/3,k=i%3;
        return backbone_atom(j).Position()[k];
    }
    i-=9;
    if (i<(3*naa)) return bb.torsional_angle(i);
    i-=(3*naa);
    if (i<numRTdof()) return RTdof(i);
    return 0;
}

void Protein::set_coord(int i, double vl)
{
    if (i<9) {
        int j=i/3,k=i%3;
        switch (k) {
          case 0: backbone_atom(j).Pos().x(vl); break;
          case 1: backbone_atom(j).Pos().y(vl); break;
          case 2: backbone_atom(j).Pos().z(vl); break;
          default: ;
        };
    } else if (i<(3*naa+9)) bb.torsional_angle(i-9,vl);
    else if (i<(3*naa+9+numRTdof())) RTdof(i-3*naa-9,vl);
}

void Protein::get_dof(vector<double> &vdof) const
{
    int ntotdof=numBBdof() +numRTdof() +9;
    vdof.resize(ntotdof,0.0);
    int i,j;

    for (int k=0;k<3;++k) {
        Vector3 v=backbone_atom(k).Position();
        vdof[3*k+0]=v.x();
        vdof[3*k+1]=v.y();
        vdof[3*k+2]=v.z();
    }

    for (i=0;i<numBBdof();++i) vdof[9+i]=BBdof(i);

    for (j=0;j<numRTdof();++j) vdof[9+i+j]=RTdof(j);
}

void Protein::get_coord(vector<double> &vdof) const
{
    int ntotdof=3*naa +numRTdof() +9;
    vdof.resize(ntotdof,0.0);
    int i,j;

    for (int k=0;k<3;++k) {
        Vector3 v=backbone_atom(k).Position();
        vdof[3*k+0]=v.x();
        vdof[3*k+1]=v.y();
        vdof[3*k+2]=v.z();
    }

    for (i=0;i<(3*naa);++i) vdof[9+i]=bb.torsional_angle(i);

    for (j=0;j<numRTdof();++j) vdof[9+i+j]=RTdof(j);
}

void Protein::set_dof(vector<double>::iterator c1, vector<double>::iterator c2)
{
    for (int k=0;k<3;++k) {
        double vx,vy,vz;
        vx=*c1++;
        vy=*c1++;
        vz=*c1++;
        backbone_atom(k).Pos(Vector3(vx,vy,vz));
    }

    for (int i=0;(i<numBBdof()) && (c1!=c2);++i,++c1) BBdof(i,*c1);

    for (int j=0;(j<numRTdof()) && (c1!=c2);++j,++c1) RTdof(j,*c1);
}

void Protein::set_coord(vector<double>::iterator c1, vector<double>::iterator c2)
{
    for (int k=0;k<3;++k) {
        double vx,vy,vz;
        vx=*c1++;
        vy=*c1++;
        vz=*c1++;
        backbone_atom(k).Pos(Vector3(vx,vy,vz));
    }

    for (int i=0;(i<(3*naa)) && (c1!=c2);++i,++c1) bb.torsional_angle(i,*c1);

    for (int j=0;(j<numRTdof()) && (c1!=c2);++j,++c1) RTdof(j,*c1);
}

bool Protein::isSameTypeAs(const Protein &p2)
{
    for (size_t i=0;i<seq.size();++i) {
        if (seq[i]!=p2.seq[i]) return false;
    }

    return true;
}

int Protein::guess_missing_coordinates(vector<bool> &specified)
{
//Stage 1: Backbone check
    bool complete=true;

    for (int ibb=0;ibb<bb.numAtoms() && complete; ++ibb) {
        complete= complete && specified[bb.atom(ibb).UniqueId()];
    }

    if (!complete) {
        if (fix_broken_backbone(specified)==0) {
            prf::cerr<<"Missing backbone coordinates could not be filled in. "
            <<"This would lead to incomplete reconstruction.\n";
            return 0;
        }
    }

    bb.reconst_bond_vectors();

//Stage 2: Backbone attachments
//H-alpha,C-beta/H-alpha2

    for (int iaa=0;iaa<naa;++iaa) {
        int ib1=0,ib2=0,ib3=0,icur=0;

        if (specified[ib1=bb.atom(3*iaa).UniqueId()] &&
            specified[ib2=bb.atom(3*iaa+1).UniqueId()] &&
            specified[ib3=bb.atom(3*iaa+2).UniqueId()]) {
            if (!(specified[icur=AA(iaa)->Hca().UniqueId()])) {
                Vector3 vcur=bb.atom(3*iaa+1).Position();
                vcur=vcur+AminoAcid::locate_Ha.fresh_eval(*bb.Bond(3*iaa+2),
                        *bb.Bond(3*iaa+1), AminoAcid::phHa);
                AtomCoordinates::vec(icur,vcur.x(),vcur.y(),vcur.z());
            }

            if (!(specified[icur=AA(iaa)->Cbeta().UniqueId()])) {
                Vector3 vcur=bb.atom(3*iaa+1).Position();

                if (AA(iaa)->OLC()!=G) {
                    vcur=vcur+AminoAcid::locate_Cb.fresh_eval(*bb.Bond(3*iaa+2),
                            *bb.Bond(3*iaa+1), AminoAcid::phCb);
                } else {
                    vcur=vcur+AminoAcid::locate_Ha2.fresh_eval(*bb.Bond(3*iaa+2),
                            *bb.Bond(3*iaa+1), AminoAcid::phHa2);
                }

                AtomCoordinates::vec(icur,vcur.x(),vcur.y(),vcur.z());
            }
        }

        //H attached to the N, if it is not a proline and not the first amino-acid
        if (iaa>0 && AA(iaa)->OLC()!=P && AA(iaa)->OLC()!=DPR &&
            specified[bb.atom(3*iaa-1).UniqueId()] &&
            specified[ib1=bb.atom(3*iaa).UniqueId()] &&
            specified[ib2=bb.atom(3*iaa+1).UniqueId()]) {
            if (!(specified[icur=AA(iaa)->labeled_atom(" H  ")->UniqueId()])) {
                Vector3 vcur=bb.atom(3*iaa).Position();
                vcur=vcur+AminoAcid::locate_Hn.fresh_eval(*bb.Bond(3*iaa),
                        *bb.Bond(3*iaa+1), pi);
                AtomCoordinates::vec(icur,vcur.x(),vcur.y(),vcur.z());
            }
        }

        //O attached to teh C, if it is not the last amino-acid
        if (iaa<(naa-1) && specified[bb.atom(3*iaa+1).UniqueId()] &&
            specified[ib1=bb.atom(3*iaa+2).UniqueId()] &&
            specified[ib2=bb.atom(3*iaa+3).UniqueId()]) {
            if (!(specified[icur=AA(iaa)->Oc().UniqueId()])) {
                Vector3 vcur=bb.atom(3*iaa+2).Position();
                vcur=vcur+AminoAcid::locate_Oc.fresh_eval(*bb.Bond(3*iaa+2),
                        *bb.Bond(3*iaa+3), pi);
                AtomCoordinates::vec(icur,vcur.x(),vcur.y(),vcur.z());
            }
        }
    }

    return 1;
}

int Protein::fix_broken_backbone(vector<bool> &specified)
{
    prf::cerr<<"Unable to repair backbone.\n";
    return 0;
}

double Protein::helix_fraction()
{
    double tot=0;

    for (int iaa=1;iaa<naa-1;++iaa) {
        if (AA(iaa)->is_helical()) ++tot;
    }

    return tot/=(naa-2);
}

double Protein::strand_fraction()
{
    double tot=0;

    for (int iaa=1;iaa<naa-1;++iaa) {
        if (AA(iaa)->is_strand()) ++tot;
    }

    return tot/=(naa-2);
}

int Protein::importXYZ(std::list<AtomRecord>::iterator istart,
                             std::list<AtomRecord>::iterator iend,
                             std::vector<bool> &assignments, int at_ligand)
{
    if (istart==iend) return 0;

    std::list<AtomRecord>::iterator lgbeg,lgend;

    lgbeg=lgend=istart;

    int curlg=lgbeg->descriptor().iresrel;

    while (lgend!=iend) {
        ++lgend;

        if (lgend==iend or lgend->descriptor().iresrel!=curlg) {
            memberLigand(at_ligand++)->read_coordinates(lgbeg,lgend, assignments);
            lgbeg=lgend;

            if (lgbeg!=iend) curlg=lgbeg->descriptor().iresrel;
        }
    }

    return 0;
}

int Protein::calc_torsions(std::vector<bool> &specified)
{
    int nerr=0;
    if ((nerr=bb.calc_torsions(specified))!=0) {
        prf::cerr<<"Failed to get "<<nerr
                <<" backbone angles from Cartesian coordinates.\n";
    }
    if (NtLigand!=NULL) {
        if (specified[NtLigand->atom(0).UniqueId()] &&
            specified[iN(0)]&&
            specified[iCa(0)]&&
            specified[iC(0)]) {
            Vector3 v1=AAv[0]->Nitrogen().Pos()-NtLigand->atom(0).Pos();
            Vector3 v2=AAv[0]->Calpha().Pos()-AAv[0]->Nitrogen().Pos();
            Vector3 v3=AAv[0]->Cprime().Pos()-AAv[0]->Calpha().Pos();
            bb.torsional_angle(0,v2.torsion(v1,v3));
        }
    } else {
        if (AAv[0]->labeled_atom("1H  ")!=NULL &&
            specified[AAv[0]->labeled_atom("1H  ")->UniqueId()] &&
            specified[iN(0)]&&
            specified[iCa(0)]&&
            specified[iC(0)]) {
            Vector3 v1=AAv[0]->Nitrogen().Pos()-AAv[0]->labeled_atom("1H  ")->Pos();
            Vector3 v2=AAv[0]->Calpha().Pos()-AAv[0]->Nitrogen().Pos();
            Vector3 v3=AAv[0]->Cprime().Pos()-AAv[0]->Calpha().Pos();
            bb.torsional_angle(0,v2.torsion(v1,v3)-pi);
        }
    }
    if (specified[iN(naa-1)]&&
        specified[iCa(naa-1)]&&
        specified[iC(naa-1)]&&
        specified[iC(naa-1)+1]) {
        Vector3 v1=AAv[naa-1]->Calpha().Pos()-AAv[naa-1]->Nitrogen().Pos();
        Vector3 v2=AAv[naa-1]->Cprime().Pos()-AAv[naa-1]->Calpha().Pos();
        Vector3 v3=AAv[naa-1]->Oc().Pos()-AAv[naa-1]->Cprime().Pos();
        bb.torsional_angle(3*naa-2,v2.torsion(v1,v3)-pi);
    }
    if (nerr<(bb.numAtoms()-3)) {
        for (size_t i=0;i<lg.size();++i) {
            lg[i]->calc_torsions(specified);
        }
    }
    return nerr;
}

prf_xml::XML_Node * Protein::make_xml_node()
{
    SeqBuild sbuild;
    prf_xml::XML_Node *xnd=new prf_xml::XML_Node("protein","");
    std::ostringstream ost;
    ost<<unid;
    xnd->set_attribute("id",ost.str());
    ost.str("");
    xnd->add_child_node("sequence",sbuild.make_string(seq,' ',' '));
    if (!(charged_CT && charged_NT)) {
        if (!charged_CT) ost<<"C\n";
        if (!charged_NT) ost<<"N\n";
        xnd->add_child_node("uncharged_ends",ost.str());
    }
    ost.precision(16);
    for (int i=0;i<3;++i) {
        Vector3 v=bb.atom(i).Position();
        ost<<v.x()<<" "<<v.y()<<" "<<v.z()<<"\n";
    }
    xnd->add_child_node("global_coordinates",ost.str());
    for (int i=0;i<numLigands();++i) {
        prf_xml::XML_Node *lgnode=memberLigand(i)->make_xml_node();
        if (lgnode!=NULL) {
            xnd->add_child_node(lgnode);
        }
    }
    return xnd;
}

void Protein::Write_XML(FILE *op)
{
    SeqBuild sbuild;
    fprintf(op,"<protein id=\"%d\">\n<sequence>\n", unid);
    fprintf(op,(sbuild.make_string(seq,' ',' ')+"\n").c_str());
    fprintf(op,"</sequence>\n");
    if (!(charged_CT && charged_NT)) {
        fprintf(op,"<uncharged_ends>\n");
        if (!charged_CT) fprintf(op,"C\n");
        if (!charged_NT) fprintf(op,"N\n");
        fprintf(op,"</uncharged_ends>\n");
    }

    fprintf(op,"<global_coordinates>\n");
    for (int i=0;i<3;++i) {
        Vector3 v=bb.atom(i).Position();
        fprintf(op,"%.16f  %.16f  %.16f\n",v.x(),v.y(),v.z());
    }
    fprintf(op,"</global_coordinates>\n");

    for (int i=0;i<numLigands();++i) {
        memberLigand(i)->Write_XML(op);
    }
    fprintf(op,"</protein>\n");
}

int Protein::metamorphose(prf_xml::XML_Node *px)
{
    if (px==NULL or px->name()!="protein") return 0;
    SeqBuild sbuild;
    std::vector<OneLetterCode> gseq;
    sbuild.parse(px->child("sequence")->value(),gseq);
    if (!good_sequence(gseq)) return 0;
    bool matched=(gseq.size()==seq.size());
    for (int i=0;matched && i<numLigands();++i) matched=(seq[i]==gseq[i]);
    if (matched) return 0;
    prf_xml::XML_Node * uce=px->child("uncharged_ends");
    if (uce) {
        std::string ucetxt=uce->value();
        charged_NT=(ucetxt.find("N")>=ucetxt.size());
        charged_CT=(ucetxt.find("C")>=ucetxt.size());
    }
    Allocate(unid,atoffset,gseq);
    return 1;
}

int Protein::assign_dof(prf_xml::XML_Node *px, int mismatch_strategy)
{
    if (px==NULL or px->name()!="protein") return 0;
    SeqBuild sbuild;
    std::vector<OneLetterCode> gseq;
    sbuild.parse(px->child("sequence")->value(),gseq);
    bool matched=(gseq.size()==seq.size());
    for (int i=0;matched && i<numLigands();++i) matched=(seq[i]==gseq[i]);
    if ((not matched) and mismatch_strategy==0) {
        prf::cerr<<"Protein::assign_dof()> Sequence mismatch while "
            <<"trying to assign coordinates from XML node. \n";
        prf::cerr<<"Sequence in XML node is : "
            <<sbuild.make_string(gseq)<<"\n";
        prf::cerr<<"The Protein object has sequence : "
            <<sbuild.make_string(seq)<<"\n";
        prf::cerr<<"Strict sequence matching has been requested.\n";
        prf::cerr<<"Therefore, nothing was assigned.\n";
        return 0;
    }
    prf_xml::XML_Node *glc=px->child("global_coordinates");
    if (glc) {
        std::string values=glc->value();
        std::istringstream ssin(values);
        double tmpq=0;
        int i=0;
        std::vector<double> gbx(9,0);
        while (ssin>>tmpq and i<9) gbx[i++]=tmpq;
        for (int j=0;j<3;++j) {
            bb.atom(j).Pos(Vector3(gbx[3*j],gbx[3*j+1],gbx[3*j+2]));
        }
    }
    // The following block is there to enable an easy conversion of
    // PROFASI text configuration files to XML configuration files
    prf_xml::XML_Node *crds=px->child("dof");
    if (crds) {
        std::string values=crds->value();
        std::istringstream ssin(values);
        double tmpq=0;
        vector<double> crd;
        while (ssin>>tmpq) crd.push_back(tmpq);
        set_dof(crd);
    }
    crds=px->child("coordinates");
    if (crds) {
        std::string values=crds->value();
        std::istringstream ssin(values);
        double tmpq=0;
        vector<double> crd;
        while (ssin>>tmpq) crd.push_back(tmpq);
        set_coord(crd);
    }

    // Where as the following block is what we expect to be used when
    // reading an XML file made by PROFASI
    for (size_t i=0;i<px->n_children();++i) {
        if (px->child(i)->name()=="group") {
            prf_xml::XML_Node *lgx=px->child(i);
            int grpindx=atoi(lgx->attribute("index").c_str());
            if (grpindx>=0 && grpindx<numLigands())
                memberLigand(grpindx)->set_coord_xml(lgx);
        }
    }
    return 1;
}
