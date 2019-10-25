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

#include "FFHandler.hh"
#include "Bias.hh"
#include "TorsionTerm.hh"
#include "HydrogenBond.hh"
#include "Hydrophobicity.hh"
#include "HydrophobicityPi.hh"
#include "ChargedSCInteraction.hh"
#include "LocExVol.hh"
#include "ExVol.hh"
#include "Extras/DistanceRestraints.hh"
#include "Extras/DihedralRestraints.hh"
#include <sstream>
#include "Extras/SBM.hh"

FFHandler::FFHandler() : HandlerBase()
{
    ff=NULL;
    p=NULL;
    ffname="FF08";
    par.option("force_field","force_field",1,"( FF08 or FF04 or FF14 )");
    par.option("skip_energy","skip",1,"(skip a term in the force field)");
    par.option("Tmix","Tx",1,"(mixing temperature * kB for DistanceRestraintsExp exponentional mixing)");
    par.option("scale_sbm","sbm",1,"(scale total energy of SBM potential)");
    par.option("aromatic_pi","pi",1,"(Set strength of FF15:HydrophobicityPi::catpi())");
    par.option("PiPiInteractionMode","PPmode",1,"(select mode)");
    par.option("SBM2mix","s2m",1,"(string representing which SBMs to mix if Tmix != 0, Example: 110 mixes the first two but not the third SBM)");
    Tmix=0;
    scaleSBM=1.0;
    arpi=1.0;
    cprep = 1.0;
    focalRes=-1;
    focalRes2=-1;
    focalResScale=1.0;
    coop = 0;
    ppmode = "1";
    s2m="";
    mixtype=1;
}

FFHandler::~FFHandler()
{
    if (ff) delete ff;
}

int FFHandler::init_ff()
{
	Logger blog;
    ForceField *newff=set_up_force_field(ffname);
    if (newff!=NULL) {
        if (ff!=NULL) delete ff;
        ff=newff;
        if (p!=NULL) ff->connect(p);
        
        
        ff->set_Tmix(Tmix);
        blog(3) << "*** init_ff:Tmix="<<Tmix<<"\n\n";
        std::cout<< "*** init_ff:Tmix="<<Tmix<<"\n\n";
        
        ff->set_SbmToMix(s2m);
        blog(3) << "*** init_ff:SbmToMix="<<s2m<<"\n\n";
        std::cout << "*** init_ff:SbmToMix="<<s2m<<"\n\n";
        
        ff->set_scaleSBM(scaleSBM);
        blog(3) << "*** init_ff:scaleSBM="<<scaleSBM<<"\n\n";
        std::cout << "*** init_ff:scaleSBM="<<scaleSBM<<"\n\n";
        
        ff->set_mixtype(mixtype);
        blog(3) << "*** init_ff:mixtype="<<mixtype<<"\n\n";
        std::cout << "*** init_ff:mixtype="<<mixtype<<"\n\n";
        
        ff->set_arpi(arpi);
        blog(3) << "*** init_ff:cationpi="<<arpi<<"\n\n";
        std::cout << "*** init_ff:cationpi="<<arpi<<"\n\n";
        
        ff->setMixable(mixable);
        std::cout<<"ff->setMixable(mixable);\n";
        
        return 1;
        
        
    } else {
        prf::cerr<<"FFHandler> Trying to set up force field "<<ffname
                <<" returned a NULL potential. \n";
        if (ff==NULL) {
            prf::cerr<<"There is no valid force field yet. Using this "
                    <<"state of the handler will likely lead to a crash.\n";
        } else {
            prf::cerr<<"The existing force field "<<ff->Name()<<" will not be changed. \n";
        }
    }
    return 0;
}

int FFHandler::set_ff(std::string ffn)
{
    ffname=ffn;
    return 1;
}



int FFHandler::parseCommand(InstructionString s)
{
	Logger blog;
    if (s.head()=="force_field") {
        if (args_check(s,1)){ 
        	set_ff(s.part(1));
        }
    } else if (s.head()=="skip" || s.head()=="skip_energy") {
    	if (args_check(s,1)) {
    		//std::cout << "skip energy:"<<s.part(1)<<"\n";
    		
    		skippedTerms = split(s.part(1),',');
			//for (int i=0;i<skippedTerms.size();i++)
    		//	std::cout << "   " << skippedTerms[i]<<"\n";
    		
    	}
    } else if (s.head()=="Tmix") {
    	if (args_check(s,1)) {
    		double t;
    		std::stringstream convert( std::string(s.part(1) ) );
    		convert>>t;
    		
    		std::cout << "found Tmix option: " << t << "\n";
    		Tmix = t;
    	}
    }else if (s.head()=="sgp" || s.head()=="scale_go_potential") {
    	if (args_check(s,1)) {
    		double n;
    		std::stringstream convert( std::string(s.part(1) ) );
    		convert>>n;
    		
    		std::cout << "found scaleGoPot option: " << n << "\n";
    		scaleSBM = n;
    	}
    } else if (s.head()=="cp" || s.head()=="cationpi") {
    	double n;
		std::stringstream convert( std::string(s.part(1) ) );
		convert>>n;
		
		std::cout << "found cationpi option: " << n << "\n";
		arpi = n;
    }else if (s.head()=="cpi" || s.head()=="cprep") {
		    	double n;
				std::stringstream convert( std::string(s.part(1) ) );
				convert>>n;
				
				std::cout << "found cprep option: " << n << "\n";
				cprep = n;
    } else if (s.head()=="residue" || s.head()=="res") {
    	int n;
		std::stringstream convert( std::string(s.part(1) ) );
		convert>>n;
		
		std::cout << "found residue option: " << n << "\n";
		blog(3) << "found residue option: " << n << "\n";
		focalRes = n-1;
    } else if (s.head()=="residue2" || s.head()=="res2") {
    	int n;
		std::stringstream convert( std::string(s.part(1) ) );
		convert>>n;
		
		std::cout << "found residue2 option: " << n << "\n";
		blog(3) << "found residue2 option: " << n << "\n";
		focalRes2 = n-1;
    } else if (s.head()=="scalepair" || s.head()=="sp") {
    	double n;
		std::stringstream convert( std::string(s.part(1) ) );
		convert>>n;
		
		std::cout << "found scale pair option: " << n << "\n";
		blog(3) << "found scale pair option: " << n << "\n";
		focalResScale = n;
    }else if (s.head()=="cooperativity" || s.head()=="coop") {
    	double n;
		std::stringstream convert( std::string(s.part(1) ) );
		convert>>n;
		
		std::cout << "found cooperativity option: " << n << "\n";
		blog(3) << "found cooperativity option: " << n << "\n";
		coop = n;
    }else if (s.head()=="PPmode" || s.head()=="PiPiInteractionMode") {
		
		ppmode = std::string(s.part(1) );
		std::cout << "found PPmode option: " << ppmode << "\n";
		blog(3) << "found PPmode option: " << ppmode << "\n";
    }else if (s.head()=="s2m" || s.head()=="SBM2mix") {
		
		s2m = std::string(s.part(1) );
		std::cout << "found s2m option: " << s2m << "\n";
		blog(3) << "found s2m option: " << s2m << "\n";
    }
    return 1;
}

void FFHandler::useEnergy(Energy *en)
{
    ff->add_external_term(en);
}

void FFHandler::skipEnergy(std::string enname)
{
    Logger blog;
    std::cout << "skipping...\n";
    if (ff->delete_term(enname)){
        blog(3)<<"FFHandler> User specified energy skip : "<<enname<<"\n";
    	std::cout <<"FFHandler> User specified energy skip : "<<enname<<"\n";
    }else{
        blog(3)<<"FFHandler> Removing term "<<enname<<" failed. The force field "
            <<"did not recognize "<<enname<<" as one of the terms.\n";
    	std::cout <<"FFHandler> Removing term "<<enname<<" failed. The force field "
        <<"did not recognize "<<enname<<" as one of the terms.\n";
    }
}

// This is a helper function. Scroll down for the real action!
bool FFHandler::get_closure(std::list<std::string> &lst,
                 std::list<std::string>::iterator st,
                 std::list<std::string>::iterator &nd)
{
    nd=st;
    ++nd;
    while (nd!=lst.end()) {
        if (*nd=="(") {
            if (!get_closure(lst,nd,nd)) return false;
            else continue;
        } else if (*nd==")") {
            ++nd;
            return true;
        }
        ++nd;
    }
    return false;
}

// This is another helper function. Ignore and scroll down.
void FFHandler::tokenize(std::string req, std::list<std::string> &tokens)
{
    tokens.clear();
    size_t i=0;
    std::string curtoken="";
    while (i<req.size()) {
        if (req[i]=='+' or req[i]=='-' or req[i]=='='
            or req[i]==':' or req[i]=='(' or req[i]==')') {
            if (!curtoken.empty()) tokens.push_back(curtoken);
            std::string separator="";
            separator+=req[i++];
            tokens.push_back(separator);
            curtoken.clear();
        } else curtoken+=req[i++];
    }
    if (!curtoken.empty()) tokens.push_back(curtoken);

    for (std::list<std::string>::iterator it=tokens.begin();
    it!=tokens.end();++it) {
        if (*it=="(") {
            std::list<std::string>::iterator jt=it,kt=it;
            bool closed=get_closure(tokens,it,jt);
            std::list<std::string> pars;
            if (closed) {
                pars.splice(pars.end(),tokens,++kt,jt);
                for (std::list<std::string>::iterator lt=pars.begin();
                lt!=pars.end();++lt) *it+=*lt;
                pars.clear();
            }
        }
    }
    for (std::list<std::string>::iterator it=tokens.begin();
    it!=tokens.end();++it) {
        std::list<std::string>::iterator jt=it;
        ++jt;
        if (jt!=tokens.end() && (*jt)[0]=='(') {
            *it+=*jt;
            tokens.erase(jt);
        }
    }
    for (std::list<std::string>::iterator it=tokens.begin();
    it!=tokens.end();++it) {
        std::list<std::string>::iterator jt=it;
        ++jt;
        if (jt!=tokens.end() && (*jt)[0]==':') {
            *it+=*jt;
            tokens.erase(jt);
            jt=it;
            ++jt;
            if (jt!=tokens.end()) {
                *it+=*jt;
                tokens.erase(jt);
            }
        }
    }
}

// This sets up a force field based on a given string. It converts
// the given string into a set of terms by first interpreting the
// string as a series of terms of this form:
// energy_term_name(parameters).
// It expands aliases such as FF04 and FF08 into their component
// terms. It interprets subtraction requests with "-" sign. It
// does not map the actual energy term names to the appropriate
// classes. Scroll further down to see how that is done.
ForceField * FFHandler::set_up_force_field(std::string req)
{
	Logger blog;
    ForceField *ff=NULL;
    if (req.empty()) req="FF08";
    std::list<std::string> tokens;
    tokenize(req,tokens);
    std::list<std::string>::iterator it=tokens.begin();
    std::list<std::string>::iterator jt=it;
    ++jt;
    std::string ffname;
    if (jt!=tokens.end() && *jt=="=") {
        ffname=*it;
        tokens.erase(it);
        tokens.erase(jt);
    } else ffname=req;
    for (it=tokens.begin();it!=tokens.end();++it) {
        //Convert FF aliases like FF08, FF04 to terms
    	if (*it=="FF15") {
            tokens.insert(it,"FF08:ExVol");
            tokens.insert(it,"FF08:LocExVol");
            tokens.insert(it,"FF08:Bias");
            tokens.insert(it,"FF08:TorsionTerm");
            tokens.insert(it,"FF08:HBMM");
            tokens.insert(it,"FF08:HBMS");
            tokens.insert(it,"FF15:HydrophobicityPi"); // NEW
            jt=it;
            --it;
            tokens.erase(jt);
    	}else if (*it=="FFx1") {
            tokens.insert(it,"FF08:ExVol");
            tokens.insert(it,"FF08:LocExVol");
            tokens.insert(it,"FF08:Bias");
            tokens.insert(it,"FF08:TorsionTerm");
            jt=it;
            --it;
            tokens.erase(jt);
        }else if (*it=="FFx2") {
        	tokens.insert(it,"FF08:LocExVol");
        	tokens.insert(it,"FF08:Bias");
            
            jt=it;
            --it;
            tokens.erase(jt);
        }else if (*it=="FFx3") {
        	tokens.insert(it,"FF08:ExVol");
        	tokens.insert(it,"FF08:LocExVol");
            
            jt=it;
            --it;
            tokens.erase(jt);
        }else if (*it=="FFx4") {
        	tokens.insert(it,"FF08:ExVol");
	        tokens.insert(it,"FF08:LocExVol");
	        tokens.insert(it,"FF08:Bias");
	        tokens.insert(it,"FF08:TorsionTerm");
	        tokens.insert(it,"FF08:HBMM");
            
            jt=it;
            --it;
            tokens.erase(jt);
        }else if (*it=="FF08") {
            tokens.insert(it,"FF08:ExVol");
            tokens.insert(it,"FF08:LocExVol");
            tokens.insert(it,"FF08:Bias");
            tokens.insert(it,"FF08:TorsionTerm");
            tokens.insert(it,"FF08:HBMM");
            tokens.insert(it,"FF08:HBMS");
            tokens.insert(it,"FF08:Hydrophobicity");
            tokens.insert(it,"FF08:ChargedSCInteraction");
            jt=it;
            --it;
            tokens.erase(jt);
        } else if (*it=="FF04") {
            tokens.insert(it,"FF04:ExVol");
            tokens.insert(it,"FF04:LocExVol");
            tokens.insert(it,"FF04:Bias");
            tokens.insert(it,"FF04:HBMM");
            tokens.insert(it,"FF04:HBMS");
            tokens.insert(it,"FF04:Hydrophobicity");
            jt=it;
            --it;
            tokens.erase(jt);
        }
    }
    //Do subtractions
    std::list<std::string> tosubtract;
    for (it=tokens.begin();it!=tokens.end();++it) {
        if (*it=="-" && it!=tokens.begin()) {
            jt=it;
            ++jt;
            if (jt!=tokens.end()) {
                tosubtract.splice(tosubtract.end(),tokens,jt);
            }
            jt=it;
            --it;
            tokens.erase(jt);
        }
    }
    for (it=tosubtract.begin();it!=tosubtract.end();++it) {
        jt=std::find(tokens.begin(),tokens.end(),*it);
        if (jt!=tokens.end()) tokens.erase(jt);
    }

    //In a new loop, convert remaining terms to FF terms
    ff=new ForceField();
    for (it=tokens.begin();it!=tokens.end();++it) {
        if ((*it)=="+" or (*it)=="-") continue;
        
        bool skipme = false;
        //std::cout << " *it:"<<*it<<"\n";
        for (int i=0; i< (int)skippedTerms.size();i++){
        	if ( ((std::string)*it).find(skippedTerms[i]) != std::string::npos){
        		skipme = true;
        	}
        }
            
        if (skipme){
        	std::cout << "skipped energy term: " << *it << "\n";
        	blog(3)<<"skipped energy term: " << *it << "\n";
        }else{
	        Energy *trm=new_energy_term(*it);
	        if (trm!=NULL) ff->add_term(trm);
	        else {
	            prf::cerr<<"FFHandler> Could not interpret "<<(*it)<<" as an energy "
	                    <<"term!\n";
	        }
        }
    }
    if (ff->n_terms()!=0) {
        ff->set_name(ffname);
        Logger(10)<<"Created interaction potential. Summary: \n"
                <<ff->summary()<<"\n";
    } else {
        prf::cerr<<"Error! Requested force-field \""<<req
                <<"\" could not be set up\n";
        delete ff;
        ff=NULL;
    }

    return ff;
}

// This function interprets a name of the form name(parameters) and
// returns a pointer to a newly created energy term when possible.
// If you want to add support for one newly created energy term on
// the commandline, this might be the only place where you need to
// add a new " else if " block.
Energy * FFHandler::new_energy_term(std::string ename)
{
    Energy *et=NULL;
    std::string eterm,pars;
    size_t parloc=ename.find_first_of('(');
    size_t parend=ename.find_last_of(')');
    parend=(parend<ename.size())?(parend-parloc-1):(ename.size()-parloc-1);

    if (parloc<ename.size()) {
        eterm=ename.substr(0,parloc);
        pars=ename.substr(parloc+1,parend);
        ename=eterm;
    }
    
    
    
    if (ename=="Bias" || ename=="FF08:Bias") et=new Bias();
    else if (ename=="TorsionTerm" || ename=="FF08:TorsionTerm")
        et=new TorsionTerm();
    else if (ename=="HBMM" || ename=="FF08:HBMM") et=new HBMM();
    else if (ename=="HBMS" || ename=="FF08:HBMS") et=new HBMS();
    else if (ename=="Hydrophobicity" || ename=="FF08:Hydrophobicity"){
        et=new Hydrophobicity();
    	if (focalRes  != -1) ((Hydrophobicity*)et)->setfocalRes(focalRes);
    	if (focalRes2 != -1) ((Hydrophobicity*)et)->setfocalRes2(focalRes2);
    	if (focalRes != -1 && focalRes2 != -1 && focalResScale != 1.0) 
    	    ((Hydrophobicity*)et)->setfocalResScale(focalResScale);
	}else if (ename=="HydrophobicityPi" || ename=="FF15:HydrophobicityPi"){
        et=new HydrophobicityPi();
        ((HydrophobicityPi*)et)->setCatPiStrength(arpi);
        ((HydrophobicityPi*)et)->setCatPiRepulsion(cprep);
        ((HydrophobicityPi*)et)->setMode(ppmode);
    	if (focalRes  != -1) ((HydrophobicityPi*)et)->setfocalRes(focalRes);
    	if (focalRes2 != -1) ((HydrophobicityPi*)et)->setfocalRes2(focalRes2);
    	if (focalRes != -1 && focalRes2 != -1 && focalResScale != 0.0) 
    	    ((HydrophobicityPi*)et)->setfocalResScale(focalResScale);
    	
    }else if (ename=="ChargedSCInteraction" || ename=="FF08:ChargedSCInteraction")
        et=new ChargedSCInteraction();
    else if (ename=="LocExVol" || ename=="FF08:LocExVol") et=new LocExVol();
    else if (ename=="ExVol" || ename=="FF08:ExVol") et=new ExVol();
    else if (ename=="FF04:Bias") {
        et=new Bias();
        et->set_pars("tune_to_ff04");
    } else if (ename=="FF04:HBMM") {
        et=new HBMM();
        et->set_pars("tune_to_ff04");
    } else if (ename=="FF04:HBMS") {
        et=new HBMS();
        et->set_pars("tune_to_ff04");
    } else if (ename=="FF04:Hydrophobicity") {
        et=new Hydrophobicity();
        et->set_pars("tune_to_ff04");
    } else if (ename=="FF04:LocExVol") {
        et=new LocExVol();
        et->set_pars("tune_to_ff04");
    } else if (ename=="FF04:ExVol") {
        et=new ExVol();
        et->set_pars("tune_to_ff04");
    } else if (ename=="Extras:ObsEnergy") {
        ObsEnergy *tmpoe=new ObsEnergy();
        oe.push_back(tmpoe);
        et=tmpoe;
    } else if (ename=="Extras:DistanceRestraints") {
        et=new DistanceRestraints();
    } else if (ename=="Extras:SBM") {
       et=new SBM();
       mixable.push_back(et);
       std::cout<<"mixable.push_back(et);\n";
       ((SBM*)et)->set_coop(coop);
       ((SBM*)et)->setScaleSBM(scaleSBM);
        
    } else if (ename=="Extras:DihedralRestraints") {
        et=new DihedralRestraints();
    }
    if (et!=NULL) et->set_pars(pars);
    return et;
}
