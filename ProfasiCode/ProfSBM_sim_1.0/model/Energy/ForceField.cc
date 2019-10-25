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

#include "ForceField.hh"
#include <sstream>
#include <typeinfo>
#include "Extras/DistanceRestraints.hh"
#include "Extras/SBM.hh"


using namespace prf;

ForceField::ForceField() {
    myname="InteractionPotential";
    etot=dele=0;
    scaleSBM = 1.0;
    hasNative = false;
    debugFF = false;
    mixtype = 1;
}

ForceField::ForceField(std::string ffname)
{
    myname=ffname;
    etot=dele=0;
    scaleSBM = 1.0;
    hasNative = false;
    debugFF = false;
    mixtype = 1;
}

ForceField::~ForceField()
{
    for (size_t i=0;i<incharge.size();++i) if (incharge[i]) delete incharge[i];
}

void ForceField::init()
{
	if (debugFF) std::cout<<"ForceField::init\n";
	//Tmix = 0;
	
	for (size_t i=0;i<e.size();++i) {
        if (e[i]!=NULL) {
        	
        	if ( e[i]->Name().find("SBM") != std::string::npos ){
        		hasNative = true;
        		((SBM*) e[i])->setDebug(debugFF);
        	}

            e[i]->init();
            e[i]->is_in_use(true);
            
        }
    }
	SBM_energies_tmp.resize(mixable.size(),9999.9999);
	minE.resize(mixable.size(),9999.9999);
	sigconsts.resize(mixable.size(),9999.9999);
	// by default mix all SBMs when mixing is enabled
	double mintmp = 0;
	double maxtmp = 0;
	double totalMin;
	
	for (size_t i = 0; i < mixable.size(); i++){
		nRestraints = ((SBM*)mixable[i])->getNrestraints();
		epsilon = - ((SBM*)mixable[i])->getPerRestraintEnergy();
		maxtmp += exp(0/Tmix);
		totalMin = epsilon*nRestraints;
		minE[i]  = totalMin;
		mintmp += exp(-totalMin/Tmix);
	}
	
	if (doMix){
		s2m.resize(mixable.size(),true);
		Mlow  = -Tmix*log(mintmp);
		Mhigh = -Tmix*log(maxtmp);
		
		if (debugFF){
			std::cout<<"Mlow="<<Mlow<<"; Mhigh="<<Mhigh<<"\n";
		}
	}else{
		s2m.resize(mixable.size(),false);
		
	}
	
	
	if (debugFF){
		for (size_t i = 0; i < mixable.size(); i++){
			std::cout<< "i"<<i<<": minE="<<minE[i]<<"; sigconst="<<sigconsts[i]<<"\n";
		}
	}

}

double ForceField::Mix(double logmix){
	return scaleSBM * logmix;
}
double ForceField::scaledMix(double logmix,double eminsum){
	return scaleSBM * (eminsum + 2*Mhigh) * (logmix - Mhigh) / Mlow ;
}


Energy * ForceField::term_called(std::string ename)
{
    Energy *ans=NULL;
    for (size_t i=0;i<e.size();++i)
        if (e[i]->Name()==ename) {
        ans=e[i];
        break;
    }
    return ans;
}

int ForceField::delete_term(std::string nm)
{
    int skipped=0;
    if (term_called(nm)!=NULL) {
        std::vector<Energy *> bkp=e;
        e.clear();
        for (size_t i=0;i<bkp.size();++i) {
            if (bkp[i]->Name()==nm) {
                delete bkp[i];
                bkp[i]=NULL;
                skipped=1;
            } else e.push_back(bkp[i]);
        }
    }
    std::cout << "success:"<<skipped<<"\n";
    return skipped;
}

std::string ForceField::summary()
{
    std::ostringstream sout;
    sout<<"Force field name "<<myname<<"\n";
    sout<<"Number of terms "<<e.size()<<"\n";
    for (size_t i=0;i<e.size();++i) {
        sout<<i<<": "<<e[i]->Name()<<"\n";
    }
    return sout.str();
}

void ForceField::print_contributions(prf::Output &op)
{
    etot=0;
    double eterm=0;
    for (size_t i=0;i<e.size();++i) {
        e[i]->refresh();
        op<<e[i]->Name()<<" = "<<(eterm=e[i]->Value())<<"\n";
        etot+=eterm;
    }
    op<<"Total = "<<etot<<"\n";
}

double ForceField::evaluate()
{
	if (debugFF) std::cout << "ForceField::evaluate\n";
	etot=dele=0;
    for (size_t i=0;i<e.size();++i) {
    	if ( e[i]->Name().find("SBM") == std::string::npos ){
    		etot+=e[i]->evaluate();
    	}
    }
    if (hasNative){
    	//etot+=NatConEvaluate();
    	SBM_E = SBM_Evaluate();
    	etot+=SBM_E;
    }
    return etot;
}

double ForceField::evaluateLambda(double lambda){
	if (debugFF) std::cout << "ForceField::evaluateLambda\n";
	etot=dele=0;
    for (size_t i=0;i<e.size();++i) {
    	if ( e[i]->Name().find("SBM") == std::string::npos ){
    		etot+=e[i]->evaluate();
    	}
    }
    if (hasNative){
    	//etot+=NatConEvaluate();
    	etot+=SBM_EnergyLambda(lambda);
    }
    return etot;
}

void ForceField::refresh()
{
	std::cout << "ForceField::refresh\n";
    etot=dele=0;
    for (size_t i=0;i<e.size();++i) {
    	if (e[i]->Name().find("SBM") == std::string::npos ){ // add dE of non-SBM terms
    		e[i]->refresh();
        	etot+=e[i]->value();
        	dele+=e[i]->deltaE();
    	}
    	
    }
    if (hasNative){
		//SBM_E = SBM_Evaluate();
		etot+=SBM_E;
		
	}
}

double ForceField::deltaE(Update *updt)
{
	if (debugFF) std::cout << "ForceField::deltaE\n";
    dele=0;
    for (size_t i=0;i<e.size();++i) {
    	if (e[i]->Name().find("SBM") == std::string::npos ){ // add dE of non-SBM terms
    		dele+=e[i]->deltaE(updt);
    	}
    }
    if (hasNative)
    	//dele += NatConDeltaE(updt);
    	SBM_E_prop = SBM_DeltaE(updt);
    	dele +=  SBM_E_prop - SBM_E;
    return dele;
}

double ForceField::SBM_EnergyLambda(double lambda){ 
	if (debugFF) std::cout<<"ForceField::SBM_Energy\n";
	double total = 0.0;
	double etmp,exptmp,Ego,logtmp,Egomixed;
	double exptmpsum = 0.0, eminsum =0;
	double exptmpsum_min=0,exptmp_min=0,nonMix_Ego_min=0,minMixedE=-9999.999;
	double nonMix_Ego =0;
	int natCounter = mixable.size();
	if (debugFF) std::cout<<"natCounter="<<natCounter<<"\n";
	
	for (size_t i=0; i<SBM_energies_tmp.size(); i++){
		
		etmp = SBM_energies_tmp[i];
		if (debugFF) std::cout<<"i:"<<i<<"; SBM_energies_tmp="<<etmp<<"\n";
		
		if (debugFF && doMix){
			std::cout<<"mix SBM "<<i<<"? "<<s2m[i]<<"\n";
		}
		
		nRestraints = ((SBM*)mixable[i])->getNrestraints();
		epsilon = - ((SBM*)mixable[i])->getPerRestraintEnergy();
    	if (doMix){
    		// do exponential mixing of native contact energies
    		
    		if (s2m[i]){
    		
	    		exptmp = exp(- (etmp)/Tmix);
	    		
	    		if (debugFF) exptmp_min = exp(minE[i]/Tmix);
	    		
	    		exptmpsum+=exptmp;
	    		eminsum += minE[i];
	    		if (debugFF) {
	    			exptmpsum_min += exptmp_min;
	    			std::cout<<e[i]->Name()<<" temp. native energy="<<etmp<<"\n";
	    		}
	    		
    		}else{
    			nonMix_Ego += etmp;
    			if (debugFF){
    				nonMix_Ego_min += minE[i];
    			}
    		}
    	}else{
		    Ego = lambda* etmp; // !!!!!!!!!!!!!!!!
    		total += Ego;
    		
    		if (debugFF) std::cout<<mixable[i]->Name()<<" Native energy="<<Ego<<"\n";
    	}
	    	
    	if (debugFF) std::cout<<"nRestraints="<<nRestraints<< "; eps="<<epsilon<<"; doMix?"<<doMix<<"; Tmix="<<Tmix<<";  goScale="<<scaleSBM<<"\n";
	}
	
	if ( doMix && natCounter > 1 ){
		// 'scaleSBM' is already applied with the function calls to 'Mix' and 'scaleMix'
		if (mixtype==1){
			Egomixed = scaledMix(-Tmix*log(exptmpsum), eminsum);
			if (debugFF) std::cout<<"eminsum="<<eminsum<<"; -Tmix*log(exptmpsum)="<<-Tmix*log(exptmpsum)<<"; exptmpsum="<<exptmpsum<<"\n";
		}else if (mixtype==2){
			Egomixed = Mix(-Tmix*log(exptmpsum));
		}
			
		if (debugFF){
	    	std::cout<<"MIXED Native energy="<<Egomixed<<"("<< (100*Egomixed/Mlow) <<"% of "<< Mlow<<")\n";
	    	std::cout<<"non-MIXED Native energy:"<< nonMix_Ego<<"\n";
	    }
	    total = Egomixed + nonMix_Ego;
	}
}

double ForceField::SBM_Energy(){ 
	if (debugFF) std::cout<<"ForceField::SBM_Energy\n";
	double total = 0.0;
	double etmp,exptmp,Ego,logtmp,Egomixed;
	double exptmpsum = 0.0, eminsum =0;
	double exptmpsum_min=0,exptmp_min=0,nonMix_Ego_min=0,minMixedE=-9999.999;
	double nonMix_Ego =0;
	int natCounter = mixable.size();
	if (debugFF) std::cout<<"natCounter="<<natCounter<<"\n";
	
	for (size_t i=0; i<SBM_energies_tmp.size(); i++){
		
		etmp = SBM_energies_tmp[i];
		if (debugFF) std::cout<<"i:"<<i<<"; SBM_energies_tmp="<<etmp<<"\n";
		
		if (debugFF && doMix){
			std::cout<<"mix SBM "<<i<<"? "<<s2m[i]<<"\n";
		}
		
		nRestraints = ((SBM*)mixable[i])->getNrestraints();
		epsilon = - ((SBM*)mixable[i])->getPerRestraintEnergy();
    	if (doMix){
    		// do exponential mixing of native contact energies
    		
    		if (s2m[i]){
    		
	    		exptmp = exp(- (etmp)/Tmix);
	    		
	    		if (debugFF) exptmp_min = exp(minE[i]/Tmix);
	    		
	    		exptmpsum+=exptmp;
	    		eminsum += minE[i];
	    		if (debugFF) {
	    			exptmpsum_min += exptmp_min;
	    			std::cout<<e[i]->Name()<<" temp. native energy="<<etmp<<"\n";
	    		}
	    		
    		}else{
    			nonMix_Ego += etmp;
    			if (debugFF){
    				nonMix_Ego_min += minE[i];
    			}
    		}
    	}else{
    		Ego = scaleSBM* etmp;
    		total += Ego;
    		
    		if (debugFF) std::cout<<mixable[i]->Name()<<" Native energy="<<Ego<<"\n";
    	}
	    	
    	if (debugFF) std::cout<<"nRestraints="<<nRestraints<< "; eps="<<epsilon<<"; doMix?"<<doMix<<"; Tmix="<<Tmix<<";  sbmScale="<<scaleSBM<<"\n";
	}
	
	if ( doMix && natCounter > 1 ){
		// 'scaleGoPot' is already applied with the function calls to 'Mix' and 'scaleMix' 
		if (mixtype==1){
			Egomixed = scaledMix(-Tmix*log(exptmpsum), eminsum);
			if (debugFF) std::cout<<"eminsum="<<eminsum<<"; -Tmix*log(exptmpsum)="<<-Tmix*log(exptmpsum)<<"; exptmpsum="<<exptmpsum<<"\n";
		}else if (mixtype==2){
			Egomixed = Mix(-Tmix*log(exptmpsum));
		}
			
		if (debugFF){
	    	std::cout<<"MIXED Native energy="<<Egomixed<<"("<< (100*Egomixed/Mlow) <<"% of "<< Mlow<<")\n";
	    	std::cout<<"non-MIXED Native energy:"<< nonMix_Ego<<"\n";
	    }
	    
	    total = Egomixed + nonMix_Ego;
	}
	
	if (debugFF && doMix){
		if (total < Mlow+nonMix_Ego_min) {
			exit(1);
		}
	}
	return total;
}
double ForceField::SBM_DeltaE(Update *updt){
	if (debugFF) std::cout<<"ForceField::SBM_DeltaE\n";
	for (size_t i=0; i<mixable.size(); i++){
		double dE = mixable.at(i)->deltaE(updt);
		if (debugFF){
			std::cout<<mixable.at(i)->Name()<<"\n";
			std::cout<<"value:" << mixable.at(i)->value() << "; delta:" << dE << "\n";
		}
		SBM_energies_tmp[i] = mixable.at(i)->value() + dE;
	}	
	return SBM_Energy();
}
double ForceField::SBM_Evaluate(){
	if (debugFF) std::cout<<"ForceField::SBM_Evaluate\n";
	for (size_t i=0; i<mixable.size(); i++){
		SBM_energies_tmp[i] = mixable.at(i)->evaluate();
	}	
	return SBM_Energy();
}
double ForceField::SBM_Minimum(){
	if (debugFF) std::cout<<"ForceField::SBM_Minimum\n";
	return Mlow;
}


double ForceField::deltaE(Update *updt, double maxde)
{
	if (debugFF) std::cout << "ForceField::deltaE2\n";
    
	dele=0;
    for (size_t i=1;i<e.size();++i) { // !NOTE! starts with i=1, i=0 is treated separately below (ExVol) because it is most expensive
    	
    	if ( e[i]->Name().find("SBM") == std::string::npos ){ // if energy term is NOT SBM?
    		if (debugFF) std::cout<<"Lund term: "<<e[i]->Name()<<" E="<<e[i]->value()<<"\n";
    		dele+=e[i]->deltaE(updt);
    	}
    	
    }
    
    if (debugFF)std::cout<<"dele before natcon: "<<dele<<"\n";
    
    if (hasNative){ // calc full SBM component energy of proposed update and add to dE
    	//dele += NatConDeltaE(updt);
    	SBM_E_prop = SBM_DeltaE(updt);
    	dele += SBM_E_prop - SBM_E;
    
	    if (debugFF){
	    	std::cout<<"dele after natcon: "<<dele<<"\n";
	    	std::cout<<"maxde before natcon: "<<maxde<<"\n";}
	    
	    maxde-=dele;
	    
	    if (debugFF) {
	    	std::cout<<"maxde after natcon: "<<maxde<<"\n";
	    	std::cout<<"dele before deltaEwithlimit: "<<dele<<"\n";}
	
	    dele+=e[0]->deltaEwithlimit(updt,maxde);
	    
	    if (debugFF) {
	    	std::cout<<"e[0]:"<<e[0]->Name()<<"\n";
	    	std::cout<<"dele after deltaEwithlimit: "<<dele<<"\n";
	    	std::cout<<"Total Energy:"<<etot<<"; total deltaE"<< dele<<"\n";}
    }else{
    	maxde-=dele;
    	dele+=e[0]->deltaEwithlimit(updt,maxde);
    }
    return dele;
}

void ForceField::accept(Update * updt)
{
	if (debugFF)std::cout<<"before ACCEPT, Etot="<<etot<<"\n";
    for (size_t i=0;i<e.size();++i) {
    	e[i]->Accept(updt);
    }
    etot+=dele;
    
    if (hasNative){
    	SBM_E = SBM_E_prop;
    	//enatcon += delenatcon;
    	//delenatcon=0;
    	SBM_E_prop = 0;
    }
    dele=0;
    
    if (debugFF)std::cout<<"after ACCEPT, Etot="<<etot<<"\n";
}

void ForceField::reject(Update *updt)
{
	if (debugFF)std::cout<<"REJECT, Etot="<<etot<<"\n";
    for (size_t i=0;i<e.size();++i) e[i]->Revert(updt);
    dele=0;
    
    if (hasNative) {
    	//delenatcon=0;
    	SBM_E_prop = 0;
    }
}


double ForceField::reset_total()
{
	if (debugFF)  std::cout<<"\nForceField::reset_total\n";
    Logger blog;
    double eerr=0,eb=0,ebp=0,et=0,eeps=1e-5,SBM_tmp=0;
    for (size_t i=0;i<e.size();++i) {

        if ( e[i]->Name().find("SBM") == std::string::npos ){ // if energy term is NOT SBM
        	ebp+=(eb=e[i]->value());
        	        
	        if ((eerr=fabs(eb-e[i]->evaluate()))>eeps) {
	            blog(3)<<"ForceField::reset_total()> "
	                    <<(e[i]->Name()).c_str()<<" : old = "<<eb
	                    <<", new = "<<e[i]->value()
	                    <<"\terror = "<<eerr<<"\n";
	            if (debugFF) std::cout << "ForceField::reset_total()> "
	            		<<(e[i]->Name()).c_str()<<" : old = "<<eb
	            		<<", new = "<<e[i]->value()
	            		<<"\terror = "<<eerr<<"\n";
	        }
        	if (debugFF) std::cout<<e[i]->Name()<<" energy="<<e[i]->value()<<"\n";
        	et+=e[i]->value();
        }
    }// end for loop
    if (hasNative){
    	SBM_tmp = SBM_Evaluate();
		if ((eerr=fabs(SBM_E-SBM_tmp))>eeps) {
	        blog(3)<<"ForceField::reset_total()> "
	                <<"NatConEvaluate"<<" : old = "<<SBM_E
	                <<", new = "<<SBM_tmp
	                <<"\terror = "<<eerr<<"\n";
	        if (debugFF) std::cout <<"ForceField::reset_total()> "
		            <<"NatConEvaluate"<<" : old = "<<SBM_E
		            <<", new = "<<SBM_tmp
		            <<"\terror = "<<eerr<<"\n";
	    }
	    
	    if (debugFF) std::cout<<"adding native contact energy:"<<SBM_tmp<<" to "<<et<<"\n";
	    SBM_E = SBM_tmp;
	    et += SBM_E;
    }
    if ( (eerr=fabs(etot-et))>eeps) {
        blog(3)<<"ForceField::reset_total()> Total E : old = "<<etot
                <<" sum of old terms = "<<ebp
                <<", new = "<<et
                <<"\terror = "<<eerr<<" -- Tmix=" << Tmix<<"\n";
        if (debugFF) std::cout<<"ForceField::reset_total()> Total E : old = "<<etot
        <<" sum of old terms = "<<ebp
        <<", new = "<<et
        <<"\terror = "<<eerr<<" -- Tmix=" << Tmix<<"\n";
        //exit(1);
    }
    etot=et;
    if (debugFF) std::cout<<"Etot = "<<etot<<"   Eerr = "<< eerr<<"\n";
    return eerr;
}

double ForceField::reset_total_silently()
{
	if (debugFF) std::cout << "ForceField::reset_total_silently\n";
    double et = 0;

    for (size_t i = 0;i<e.size();++i) {
    	if ( e[i]->Name().find("SBM") == std::string::npos ){
    		et += e[i]->evaluate();
    	}
    }
    if (hasNative){
    	SBM_E = SBM_Evaluate();
    	et+=SBM_E;
    }
    double eerr=fabs(etot-et);
    etot=et;
    if (debugFF){
    	std::cout<<"Eerr:"<<eerr<<"\n";
    	if (eerr>=1) exit(1);
    }
    return eerr;
}

Etot::Etot() : myff(NULL)
{
    Name("Etot");
    grdtyp=1;
}

Etot::~Etot() {}

void Etot::rangeEstimate(double &x1, double &x2)
{
    if (myff!=NULL) {
    	
        double lx1(0),lx2(0);
        x1=x2=0;
        for (size_t i=0;i<myff->n_terms();++i) {
            
        	myff->term(i)->rangeEstimate(lx1,lx2);
        	if ( myff->term(i)->Name().find("SBM") != std::string::npos ){
        		lx1 *= myff->get_scaleSBM();
        	}
            std::cout << "range estimate " <<  myff->term(i)->Name() << " " << lx1 << " " << lx2 << "\n";
            x1+=lx1;
            x2+=lx2;
        }
        lx1=0.5*(x1+x2);
        lx2=0.5*(x2-x1);
        x1=lx1-0.5*lx2;
        x2=lx1+0.5*lx2;
        // The above symmetric shrinking of range is because the different
        // energy terms normally do not achieve their extreme values
        // for the same configurations.

        if (!userbinsz) xbin0=0.5;
    }
}

double Etot::evaluate()
{
	//std::cout << "Etot::evaluate\n";
    return myff->evaluate();
}

double Etot::delta(Update *u)
{
	//std::cout << "Etot::delta\n";
    return myff->deltaE(u);
}

void Etot::refresh()
{
	//std::cout << "Etot::refresh\n";
    obsval=myff->value();
}


Ego::Ego() : myff(NULL)
{
    Name("Ego");
    grdtyp=1;
}

Ego::~Ego() {}

void Ego::rangeEstimate(double &x1, double &x2)
{
    if (myff!=NULL) {
    	
        x1=myff->SBM_Minimum();
        x2=0;
        std::cout<<"Ego range:"<<x1<<","<<x2<<"\n";
    }
}

double Ego::evaluate()
{
	//std::cout << "Etot::evaluate\n";
    return myff->SBM_Evaluate();
}

double Ego::delta(Update *u)
{
	//std::cout << "Etot::delta\n";
    return myff->SBM_DeltaE(u);
}

void Ego::refresh()
{
	//std::cout << "Etot::refresh\n";
    obsval=myff->getNatEn();
}

