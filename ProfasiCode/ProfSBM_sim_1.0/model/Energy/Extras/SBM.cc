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

#include <fstream>
#include <deque>
#include "../../Aux/fileutils.hh"
#include "SBM.hh"

using namespace prf;

SBMrestraint::SBMrestraint() : atom1(-1),atom2(-1),f(NULL),res1(0),res2(0) {
	perRestraintEnergy=0;
}

// f is deliberately not destroyed in the destructor here!
SBMrestraint::~SBMrestraint() {}

SBMrestraint::SBMrestraint(const SBMrestraint &dr)
{
	perRestraintEnergy=0;
    atom1=dr.atom1;
    atom2=dr.atom2;
    res1=dr.res1;
    res2=dr.res2;
    f=dr.f;
}

SBMrestraint & SBMrestraint::operator =(const SBMrestraint &dr)
                                                 {
    if (this!=&dr) {
        atom1=dr.atom1;
        atom2=dr.atom2;
        res1=dr.res1;
        res2=dr.res2;
        f=dr.f;
    }
    return *this;
}

void SBMrestraint::delete_restraint_function()
{
    if (f) delete f;
    f=NULL;
}

double SBMrestraint::evaluate()
{
    double dist=AtomCoordinates::dist(atom1,atom2);
    return (*f)(dist);
}

double SBMrestraint::estimate_upper_bound()
{
    double dist=0.5*AtomCoordinates::boxL();
    return (*f)(dist);
}

double SBMrestraint::estimate_lower_bound()
{
    //double dist=0.5*AtomCoordinates::boxL();
    return f->estimate_min();
}

int SBMrestraint::set_pars(prf_xml::XML_Node *pars,Population *p)
{
    prf::Logger blog(10);
    atom1=atom2=-1;
    res1=res2=-1;
    if (f!=NULL) delete f;
    f=NULL;
    if (pars==NULL or pars->name()!="restraint") return 0;
    
    if (pars->child("atom1")!=NULL) 
    	atom1=get_aid(pars->child("atom1")->value(),p,res1);
    else if (pars->child("atom_1")!=NULL)
    	atom1=get_aid(pars->child("atom_1")->value(),p,res1);
    
    
    
    if (pars->child("atom2")!=NULL) 
    	atom2=get_aid(pars->child("atom2")->value(),p,res2);
    else if (pars->child("atom_2")!=NULL)
        atom2=get_aid(pars->child("atom_2")->value(),p,res2);
    
    //
    //std::cout << "atom1="<<atom1<<", atom2="<<atom2<<"\n";
    //std::cout << "res1="<<res1<<", res2="<<res2<<"\n";
    //std::cout << "res1="<<getRes1()<<", res2="<<getRes2()<<"\n";
    
    
    std::string rftype=pars->attribute("type");
    if (rftype.empty() or rftype=="0" or rftype=="quadratic") {
        f=new RestraintFunction();
    } else if (rftype=="1" or rftype=="power_law" or rftype=="powerlaw") {
        f=new PowerLawRestraint();
    } else if (rftype=="2" or rftype=="flattened_power_law") {
        f=new FlattenedPL();
    } else if (rftype=="3" or rftype=="gaussian" or rftype=="Gaussian") {
        f=new GaussianRestraint();
    } else if (rftype=="4" or rftype=="LJ" or rftype=="lj") {
    	f=new LennardJonesNativeRestraint();
    } else if (rftype=="5" or rftype=="GAUSS" or rftype=="gauss") {
    	f=new GaussianNativeRestraint();
    } else if (rftype=="6" or rftype=="DUALGAUSS" or rftype=="dualgauss") {
    	f=new DualGaussianNativeRestraint();
    } else if (rftype=="7" or rftype=="FDUALGAUSS" or rftype=="fdualgauss") {
    	f=new FixedDepthDualGaussianNativeRestraint();
    } else if (rftype=="8" or rftype=="FMULTIGAUSS" or rftype=="fmultigauss") {
    	f=new FixedDepthMultiGaussianNativeRestraint();
    } else if (rftype=="9" or rftype=="FMULTIGSMOOTH" or rftype=="fmultigsmooth") {
    	f=new FixedDepthMultiGaussianSmoothNativeRestraint();
    	//prf::cout<<"Hooray!\n";
    }
    
    //prf::cout << " atom1:" << atom1 << " " << pars->child("atom_1")->value() << " atom2:" << atom2 << " " << pars->child("atom_2")->value() << " " << rftype << "\n";
    
    if (atom1>0 && atom2>0 && f!=NULL) {
        if (f->set_pars(pars->child("parameters"))==0) {
            prf::cerr<<"Error while trying to set parameters for "<<rftype
                    <<" between atoms "<<atom1<<" and "<<atom2<<"\n";
            //prf::cout<<"Error while trying to set parameters for "<<rftype
            //<<" between atoms "<<atom1<<" and "<<atom2<<"\n";
        } else {
            blog<<"Created "<<rftype<<" restraint between atoms with unique ids "
                    <<atom1<<" and "<<atom2<<"\n";
            //prf::cout<<"Created "<<rftype<<" restraint between atoms with unique ids "
            //<<atom1<<" and "<<atom2<<"\n";
            return 1;
        }
    }
    if (f!=NULL) delete f;
    
    
    return 0;
}

int SBMrestraint::get_aid(std::string atm,Population *p, int &r)
{
	prf::Logger blog(3);
    std::deque<std::string> parts;
    prf_utils::split_str(atm,'/',parts,4);
    int chid,resid;
    chid=atoi(parts[0].c_str());
    resid=atoi(parts[1].c_str());
    r=resid;
    if (chid<0 or chid>=p->NumberOfChains()) {
        prf::cerr<<"Invalid chain id "<<chid<<" in restraint specification"
                <<" :"<<atm<<"\n";
        return -1;
    }
    if (resid<0 or resid>=p->Chain(chid)->numLigands()) {
        prf::cerr<<"Invalid residue number "<<resid<<" in restraint "
                <<"specification :"<<atm<<"\n"
                <<"Chain "<<chid<<" has "<<p->Chain(chid)->numLigands()
                <<" groups.\n";
        return -1;
    }
    //NEW: check atom name first to see if it is _CA_
    std::string atmname=parts[3];
        
    
        
    // NEW: if atom is _CA_, then do not raise error when residue type mismatch
    std::string aa=p->Chain(chid)->memberLigand(resid)->TLC();
    if (aa!=parts[2]) {
    	if (atmname=="_CA_"){
    		blog << "Group "<<resid<<" in chain "<<chid<<" is a "<<aa
            <<" where as it is given as a "<<parts[2]<<" in restraint "
            <<"specification :"<<atm<<" -- IGNORING SINCE ATOM IS _CA_ !!!\n";
    	}else{
	        prf::cerr<<"Group "<<resid<<" in chain "<<chid<<" is a "<<aa
	                <<" where as it is given as a "<<parts[2]<<" in restraint "
	                <<"specification :"<<atm<<"\n";
	        return -1;
    	}
    }
    
    for (size_t i=0;i<atmname.size();++i) if (atmname[i]=='_') atmname[i]=' ';
    Atom *tmpatm=p->Chain(chid)->memberLigand(resid)->labeled_atom(atmname);
    if (tmpatm!=NULL) {
    	//std::cout << "found" << tmpatm->UniqueId() << "\n";
        return tmpatm->UniqueId();
    }
    //std::cout << "failed" << tmpatm->UniqueId() << "\n";
    return -1;
    
}

SBM::SBM() : Energy(), nRestraints(0), nchanges(0)
{
    Name("SBM");
    normEnergy = 0;
    perRestraintEnergy = 0;
    coop = 1.0;
    debug = false;
    scale_SBM = 1;
}
SBM::SBM(std::string m) : Energy(), nRestraints(0), nchanges(0)
{
    Name("SBM");
    normEnergy = 0;
    perRestraintEnergy = 0;
    coop = 1.0;
    posmapfilename=m;
    debug = false;
    scale_SBM = 1;
}


SBM::SBM(double ne) : Energy(), perRestraintEnergy(0), nRestraints(0), nchanges(0)
{
    Name("SBM");
    normEnergy = ne;
    coop = 1.0;
    debug = false;
    scale_SBM = 1;
}

SBM::~SBM()
{
    for (size_t i=0;i<c.size();++i) {
        c[i].delete_restraint_function();
    }
}

void SBM::set_pars(std::string s)
{
    filename=s;
}



void SBM::init()
{
    prf_xml::XML_Node *root=prf_xml::get_xml_tree(filename);
    if (root==NULL or root->name()!="sbm") {
        if (root==NULL) {
            prf::cerr<<Name()<<"> No valid XML tree in "<<filename<<"\n";
        } else {
            prf::cerr<<Name()<<"> Root node in "<<filename
                    <<" is not <sbm>\n";
            delete root;
        }
        return;
    }
    
    restraintID="";
    
    typedef std::map<std::string,std::string>::iterator it;
    for(it iter = root->attributes().begin(); iter != root->attributes().end(); ++iter) {
    	if (iter->first == "rid")
    		restraintID = std::string(root->attribute("rid"));
    		std::string name = "SBM"+ std::string(restraintID);
    		Name(name);
        // iterator->first = key
        // iterator->second = value
        // Repeat if you also want to iterate through the second map.
    }
    
    
    root->interpret_formatted_data();
    std::deque<SBMrestraint> tmp;
    prf::Logger blog(3);
    
    for (size_t i=0;i<root->n_children();++i) {
        prf_xml::XML_Node *rst=root->child(i);
        if (rst->name()!="restraint") continue;
        SBMrestraint res;
        if (res.set_pars(rst,p)==1){
        	tmp.push_back(res);
        	//prf::cout <<"getRes... "<< res.getRes1() <<" " << res.getRes2()<< " "<< "\n";
        	//prf::cout << tmp.size() << "\n";
        }
    }
    c.resize(tmp.size());
    
    nRestraints = tmp.size();
    if (false){
    	for (int i=0; i < nRestraints;i++){
    		if (i>0 && perRestraintEnergy != scale_SBM * ((SBMrestraint) tmp[i]).getPerRestraintEnergy()){
    			std::cout<<"ERROR: epsilons of SBM have inhomogeneous values!\n";
    			exit(1);
    		}
    		perRestraintEnergy = scale_SBM * ((SBMrestraint) tmp[i]).getPerRestraintEnergy();
    		((SBMrestraint) tmp[i]).setPerRestraintEnergy(perRestraintEnergy);
    	}
    }else{
    	perRestraintEnergy = ((SBMrestraint) tmp[0]).getPerRestraintEnergy();
    }
    std::cout<<"SBMrestraint::perRestraintEnergy = "<<perRestraintEnergy<<"\n";
    c.assign(tmp.begin(),tmp.end());
    //prf::cout << c.size() << "\n";
    //for (int i=0; i<c.size();i++){
    //	prf::cout <<"what's in c? " << ((NativeRestraint)c[i]).getRes1() <<" " << ((NativeRestraint)c[i]).getRes2()<<"\n";
    //	
   //}
    blog<<Name()<<"> Created "<<c.size()<<" SBM restraints.\n";
    if (root) delete root;
}

double SBM::evaluate()
{
    vval=delv=0;
    for (size_t i=0;i<c.size();++i) {
        vval+=c[i].evaluate();
    }
    //if (coop != 1.0 && vval < 0.0){
    //	vval = - fabs( pow(-vval,coop)/pow(nRestraints,coop-1));
    //}
    return vval;
}

double SBM::deltaE(Update *updt)
{	
	// assuming only backbone atoms restrained!!
	if (debug) std::cout<<"SBM::deltaE(Update *updt) - "<<updt->Name()<<"\n";
			
	if (updt->sidechain_update()==true){
		delv=0;
		
		
	} else if (updt->Name()=="Pivot"){// || updt->Name()=="Rotation" || updt->Name()=="Translation"){
		delv=0;
		double total=0;
		double delv1=0;
		double delv1rest=0;
		int s = updt->n_residue_rigid_ranges();
		//prf::cout << "n rng:"<<s<<"\n";
		std::vector<std::pair<int,int> > *rng = updt->residue_rigid_ranges();
		for (size_t i=0;i<c.size();++i) {
			//prf::cout << "res1,2:"<<c[i].getRes1()<<","<<c[i].getRes2()<<"\n";
			bool isrig;
			int first,last =-1;
				
				
			for (int it=0; it<s;it++) {
				
				updt->residue_rigid_range( it,  first, last);
				
				//QUICK HACK FOR PIVOT UPDATES, ONLY WORKS FOR MONOMERS!
				// IS THIS A BUG IN PROFASI? WHY IS ONE RIGID RANGE ALWAYS ONE RESIDUE LONG IN PIVOT???
				// ONE-RESIDUE RIGID ELEMENT IS THE PIVOT ELEMENT - RANGES BELOW AND ABOVE ARE RIGID

				bool isrig = false;
				if (first==last){
					if (it==0){
						first = p->chain_start(0);
						last = last-1;
					}else if (it==1){
						last = p->chain_end(0)-1;
						first = first+1;
					}
				}
					
				//prf::cout <<s<< " " <<i<<" "<< c[i].getRes1()<< " "<< c[i].getRes2()<< " "<< it << " " << first << " " << last << "\n";
				
				if (c[i].is_in_range(first, last)){
					//prf::cout <<"true\n";
					isrig = true;
				}
			}
			
			
			if (! isrig)//isContactInRigidRange(i,s, rng, updt))
				delv1 += c[i].evaluate();
			//else
			//	delv1rest += c[i].evaluate();
			
			//total += c[i].evaluate();
		}
		double eold=vval;
		
		//double delv2= evaluate();
		//double delv2=(vval-eold);
		vval=eold;
		delv=(delv1-vval);
		//if (delv != 0.0)
		//	prf::cout << "   Exhaustive calculation: vval (total contact energy before update):"<<vval<<"\n";
		/*
		if (vval != delv2){
			prf::cout << updt->Name() << " side chain:" << updt->sidechain_update()<<"\n";
			prf::cout << "   Exhaustive calculation: vval (total contact energy before update):"<<vval<<" exact calculation of new energy :"<< delv2 <<"\n";
			prf::cout << "   Fast calculation:       delv1 (new energy of altered chain)="<<delv1<<" delv (energy change after proposed update):"<<delv << "(unaltered chain energy:"<<delv1rest<<", sum:"<<total <<") \n";
		}*/
	} else {// This is the trivial way to calculate delta. Normally, there
	    // is a better way to calculate it based on the properties of the
	    // particular energy term and the conformational move represented
	    // by the update updt.
		double eold=vval;
	    evaluate();
	    delv=(vval-eold);
	    vval=eold;
	    
	    
	    	//prf::cout<< "delv="<<delv<<"\n";
	    //prf::cout << updt->Name() << " side chain:" << updt->sidechain_update() << "\n";
	}
	
	
	
	return delv;
}

void SBM::rangeEstimate(double &x1,double &x2)
{
    x1=0;
    x2=0;
    for (size_t i=0;i<c.size();++i) {
        x1+=c[0].estimate_lower_bound();
    }
}

