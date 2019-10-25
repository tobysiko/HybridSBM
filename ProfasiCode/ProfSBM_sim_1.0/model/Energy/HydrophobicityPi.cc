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
#include <stdlib.h>
#include <cstdlib>
#include "HydrophobicityPi.hh"
#include <iostream>
#include <sstream>
#include <math.h>
#include <cassert>
#include <algorithm>

using namespace std;

namespace prf
{

	const std::string HydrophobicityPi::pipimap_filenames[] = { "FF.pidat",
																"YY.pidat",
																"WW.pidat",
																"FY.pidat",
																"FW.pidat",
																"YW.pidat",
																"FR.pidat",
																"YR.pidat",
																"WR.pidat",
																"FK.pidat",
																"YK.pidat",
																"WK.pidat"}; // files containing geometric pi-pi parameters
    const int HydrophobicityPi::nPiFilenames = 12;
	const int HydrophobicityPi::nHPgrp=11;
    const double HydrophobicityPi::hpstr[]={ 0.3, 0.4, 0.4, 0.6, 0.8, 0.8, 0.8, 0.8, 1.1, 1.6, 1.6 };

    string HydrophobicityPi::hatms[][6]={ // for contact-based hydrophobicity
        {" CB "," CG ","    ","    ","    ","    "}, //R
        {" CB "," CG "," SD "," CE ",    "",    ""}, //M
        {" CB "," CG "," CD ","    ","    ","    "}, //K
        {" CB "," CG1"," CG2",    "",    "",    ""}, //V
        {" CB "," CG1"," CG2"," CD1",    "",    ""}, //I
        {" CB "," CG "," CD1"," CD2",    "",    ""}, //L
        {" CB "," CG "," CD ",    "",    "",    ""}, //P
        {" CB "," CG "," CD ",    "",    "",    ""}, //DPR
        {" CG "," CD1"," CE1"," CZ "," CE2"," CD2"}, //Y
        {" CG "," CD1"," CE1"," CZ "," CE2"," CD2"}, //F
        {" CG "," CD1"," CD2"," CE3"," CZ3"," CH2"}  //W
    };
    
    int HydrophobicityPi::htyp(OneLetterCode cod)
    {
        OneLetterCode hpcdord[]={NONE,R,M,K,V,I,L,P,DPR,Y,F,W};
        int i=0;

        for (i=nHPgrp;i>0&&hpcdord[i]!=cod;--i);

        return i-1;
    }
    
    
    // PI-PI interaction
    const int HydrophobicityPi::nPPgrp=3;
    const double HydrophobicityPi::ppstr[]={ 1.0, 1.0, 1.0 };
    
    
    int HydrophobicityPi::ptyp(OneLetterCode cod)
    {
        OneLetterCode ppcdord[]={NONE,Y,F,W};
        int i=0;

        for (i=nPPgrp;i>0&&ppcdord[i]!=cod;--i);

        return i-1;
    }
    
    
    ///////
    
    string HydrophobicityPi::cenatms[][3]={ // for aromatic pi interactions
            {" CG "," CE1"," CE2"}, //Y
            {" CG "," CE1"," CE2"}, //F
            {" CD2"," CZ3"," CZ2"}  //W
            //{" CD1"," CZ2"," CE3"}  //W
        };
    
    const int HydrophobicityPi::nCPgrp=2;
    const double HydrophobicityPi::cpstr[]={ 1.0, 1.0, 1.0 };
    
    int HydrophobicityPi::ctyp(OneLetterCode cod)
    {
        OneLetterCode ppcdord[]={NONE,K,R};
        int i=0;

        for (i=nCPgrp;i>0&&ppcdord[i]!=cod;--i);

        return i-1;
    }
    
    
    ///////
    string HydrophobicityPi::catatms[][2]={ // for cation-pi interactions
        {" NZ "," CE "}, //K
        {" CZ "," NE "}  //R
    };
    
    
    double vectorNorm(Vector3 v){
        	return sqrt( abs(v.x())*abs(v.x()) + abs(v.y())*abs(v.y()) + abs(v.z())*abs(v.z()) );
        	
        }
    
    Vector3 HydrophobicityPi::centroid(Ligand* lig){
    	
    	OneLetterCode tmp_olc = lig->OLC();
    	int ty = ptyp(tmp_olc);
    	
    	assert(ty!=-1);
    	
    	string* centroid_atoms = cenatms[ty];
    	
    	Vector3 com3 = AtomCoordinates::CenterOfMass3(lig->labeled_atom(centroid_atoms[0])->UniqueId(), 
				  lig->labeled_atom(centroid_atoms[1])->UniqueId(), 
				  lig->labeled_atom(centroid_atoms[2])->UniqueId());
    	return com3;
    }
    
    
    Vector3 HydrophobicityPi::displacePoint(Vector3 point, Vector3 vect, double distance){
    	double vlen = vectorNorm(vect);
    	Vector3 unitV = vect * (1.0/vlen);
    	Vector3 newVect = unitV * distance;
    	Vector3 newPoint = point + newVect;
    	return newPoint;
    }
    
    Vector3 HydrophobicityPi::normal(Ligand* lig){
    	OneLetterCode tmp_olc = lig->OLC();
    	int ty = ptyp(tmp_olc);
    	
    	assert(ty!=-1);
    	
    	string* planar_atoms = cenatms[ty];
    	
    	Vector3 v1 = AtomCoordinates::vec(lig->labeled_atom(planar_atoms[0])->UniqueId());
    	
    	Vector3 v2 = AtomCoordinates::vec(lig->labeled_atom(planar_atoms[1])->UniqueId());
    	Vector3 v3 = AtomCoordinates::vec(lig->labeled_atom(planar_atoms[2])->UniqueId());
    	
    	Vector3 v1v2 = v2 - v1;
    	Vector3 v1v3 = v3 - v1;
    	
    	Vector3 vcross = v1v2.cross(v1v3);
    	
    	return vcross;
    }
    
    HydrophobicityPi::HydrophobicityPi() : Energy()
    {
        Name("HydrophobicityPi");
        a2=3.7*3.7;
        b2=4.5*4.5;
        d2=b2-a2;
        slpe=-1.0/d2;
        strnn=0.0;
        strnnn=0.5;
        sclfct=1.0;
        grdtyp=2;
        
        rmin = 3.0;
        rmax = 12.0;
        thetamin = 0.0;
        thetamax = 90.0;
        phimin = 0.0;
        phimax = 90.0;
        rstep = 0.1; // 90 bins
        thetastep = 1.0; // 90 bins
        phistep = 1.0; // 90 bins
        
        arpistrength = 1.0;
        
        debug = false;
        
        mode = "1";
    }

    HydrophobicityPi::~HydrophobicityPi() {}

    double HydrophobicityPi::InterChain(int ich, int jch)
    {
        double e=0;
        int hi,hj,jbeg;

        for (int i=p->chain_start(ich);i<p->chain_end(ich);++i) {
            if ((hi=ih[i])<0) continue;

            if (ich==jch) jbeg=i+1;
            else jbeg=p->chain_start(jch);

            for (int j=jbeg;j<p->chain_end(jch);++j) {
                if ((hj=ih[j])<0) continue;

                e+=Mehp[NHAA*hi+hj];
            }
        }

        return sclfct*e;
    }

    double HydrophobicityPi::interhp()
    {
        double e=0;

        for (int ich=0;ich<NC();++ich) {
            for (int jch=ich+1;jch<NC();++jch) {
                e+=InterChain(ich,jch);
            }
        }

        return e;
    }

    void HydrophobicityPi::init()
    {
        if (initialized) return;
        Logger blog;
        double tmpx;
        N=p->NumberOfLigands(),NHAA=0;

        for (int i=0;i<p->NumberOfLigands();++i) {
            if (htyp(p->ligand(i)->OLC())>=0) ++NHAA;
        }

        blog(15)<<Name()<<"> Counted "<<NHAA<<" hydrophobic ligands\n";

        //NHAA = total number of  hydrophobic amino acids or endgroups
        //In the following, the strange initialisations are an attempt to
        //try and crash the program if something illegal happens.
        strength.resize(NHAA*NHAA,0.0);
        atom.resize(NHAA*6,-200000);
        natom.resize(NHAA,0);
        soften.resize(NHAA,0);
        ih.resize(N,-1);
        seq.resize(N,NONE);
        hlig.resize(NHAA,-1);
        chainof.resize(NHAA,-1);
        Mehp.resize(NHAA*NHAA,0);
        Vehp.resize(NHAA*NHAA,0);
        changed.resize(NHAA*NHAA,0);
        centroids.resize(NHAA,NONE);
        normals.resize(NHAA, NONE);
        holc.resize(NHAA,NONE);
        ispi.resize(NHAA,NONE);//aromatics /pi
        iscat.resize(NHAA,NONE);//cations
        Ligand *lg;
        
        initPiInteractionMap();
        
        //int ligind = 0;
        
        prf::cout << "HydrophobicityPi::init\n";
        
        for (int i1=0,k1=0,l1=0;i1<NC();++i1) {

            for (int j1=0;j1<p->Chain(i1)->numLigands();++j1,++k1) {
                lg=p->Chain(i1)->memberLigand(j1);

                if (htyp(seq[k1]=lg->OLC()) ==-1) continue;
                hlig[l1] = lg->UniqueId(); // for each hydrophobic, store ligand index
                holc[l1] = lg->OLC();
                
                if (lg->OLC()==Y || lg->OLC()==W || lg->OLC()==F){
                	ispi[l1] = true;
                	iscat[l1] = false;
                }else if (lg->OLC()==K || lg->OLC()==R){
                	iscat[l1] = true;
                	ispi[l1] = false;
                
                }else{
                	ispi[l1] = false;
                	iscat[l1] = false;
                }
                
                
                prf::cout << "HydrophobicityPy::init "<< Groups::mapOLC2Char(lg->OLC()) << "," << lg->UniqueId()+1 << " - "<<l1<< " Pi?"<<ispi[l1]<< " Cat?"<< iscat[l1] <<"\n";
                
                
                for (int i2=0,k2=0,l2=0;i2<NC();++i2) { //chains
                    for (int j2=0;j2<p->Chain(i2)->numLigands();++j2,++k2) { //ligands per chain
                        if (htyp(seq[k2]=p->Chain(i2)->memberLigand(j2)->OLC())==-1) continue; // check if OLC in htyp

                        tmpx=hpstr[htyp(seq[k1])] + hpstr[htyp(p->ligand(k2)->OLC())]; // interaction strength btw ligands k1 and k2

                        if (i1==i2&&(k2==(k1+1)||k2==(k1-1))) tmpx*=strnn; // nearest neighbours have 0.0 interaction

                        if (i1==i2&&(k2==(k1+2)||k2==(k1-2))) tmpx*=strnnn; // next nearest have 0.5 interaction

                        strength[l1*NHAA+l2++]=tmpx;
                    }
                }

                chainof[l1]=i1;

                if (seq[k1]==P || seq[k1]==F || seq[k1]==Y
                    || seq[k1]==W || seq[k1]==DPR) soften[l1]=1;
                else if (seq[k1]==L || seq[k1]==I || seq[k1]==V) soften[l1]=2;

                ih[k1]=l1++;
            }
        }

        blog(50)<<"Hydrophobic Ligands list ...\n";

        for (int i=0;i<p->NumberOfLigands();++i) {
            if (ih[i]>=0) {
                blog<<ih[i]<<"  "<<Groups::mapOLC2Char(seq[i])<<" on chain "
                <<chainof[ih[i]]<<" index "<<p->ligand(i)->SeqSerial()
                <<"\n";
            }
        }

        blog<<"Pairwise hydrophobicity strength matrix ...\n";

        for (int i=0;i<NHAA;++i) {
            for (int j=0;j<NHAA;++j) {
                blog<<strength[i*NHAA+j]<<"  ";
            }

            blog<<"\n";
        }

        for (int ilg=0,n=0;ilg<p->NumberOfLigands();++ilg) {
            
        	if ((n=ih[ilg]) < 0) continue; // if ligand not hydrophobic, skip

            for (int i=0;i<6;++i) {
                Atom *at1=p->ligand(ilg)->labeled_atom(hatms[htyp(seq[ilg])][i]);

                if (at1!=NULL) {
                	atom[6*n+i]=at1->UniqueId(); // write all atoms of all ligands into single linear array
                	natom[n]++; // number of atoms per ligand
                }
            }
        }

        double max_at_dist=(14+sqrt(a2));

        max_at_dist2=max_at_dist*max_at_dist;
        max_at_dist=AtomCoordinates::boxL()-max_at_dist;
        max_cmp_dst2=max_at_dist*max_at_dist;
        
        blog(3)<<"PiPi Interaction Mode:"<<mode<<"\n";
        std::cout<<"PiPi Interaction Mode:"<<mode<<"\n";
        
        initialized=true;
    }

    void HydrophobicityPi::DisplayMatrix()
    {
        for (int i=0;i<NHAA;++i) {
            for (int j=0;j<NHAA;++j) prf::cout <<Mehp[NHAA*i+j]<<"\t";

            prf::cout <<"\n";
        }

        for (int i=0;i<NHAA;++i) {
            for (int j=0;j<NHAA;++j) prf::cout <<Vehp[NHAA*i+j]<<"\t";

            prf::cout <<"\n";
        }
    }

    double HydrophobicityPi::evaluate()
    {
        double eterm=0;
        vval=delv=0;
        nchanges=0;

        for (int i=0;i<NHAA;++i) {
            for (int j=i+1;j<NHAA;++j) {
                eterm=hp_pair(i,j,chainof[i]==chainof[j]?
                              AtomCoordinates::s2:
                              AtomCoordinates::dist2);
                vval+=eterm;
                Mehp[NHAA*i+j]=Mehp[NHAA*j+i]=eterm;
            }
        }

        return vval*=sclfct;
    }

    double HydrophobicityPi::gradientXYZ(std::valarray<double> &gx)
    {
        double eterm=0;
        vval=delv=0;
        nchanges=0;
        gx=0;
        for (int i=0;i<NHAA;++i) {
            for (int j=i+1;j<NHAA;++j) {
                eterm=dhp_pair(i,j,chainof[i]==chainof[j],gx);
                vval+=eterm;
                Mehp[NHAA*i+j]=Mehp[NHAA*j+i]=eterm;
            }
        }

        return vval*=sclfct;

    }

    void HydrophobicityPi::rangeEstimate(double &x1, double &x2)
    {

        x1=-1.2*sclfct*max(NHAA,1);

        x2=0;
        //mere guess.
    }

    double HydrophobicityPi::deltaE(Update *updt)
    {
        delv=0;
        nchanges=0;

        if (updt->sidechain_update()) {
            int ind=ih[updt->change(0).info.group];

            if (ind<0 || atom[6*ind+natom[ind]-1]<updt->begin_atom() ||
                atom[6*ind]>=updt->end_atom()) return 0;
        }

        double eterm=0,de=0;

        updt->n_residue_rigid_ranges();
        vector<pair<int,int> > *rng=updt->residue_rigid_ranges();

        for (int it=0; it<updt->n_residue_rigid_ranges();++it) {
            int ires1=(*rng)[it].first, ires2=(*rng)[it].second+1;

            for (int jres=0,prng=0;jres<N;) {
                int ihres=0,jhres=0;

                if ((jhres=ih[jres])<0) {++jres;continue;} // residue j is not hydrophobic, skip.

                if (prng<=it && jres>=(*rng)[prng].first) {
                    jres=(*rng)[prng++].second+1;
                    continue;
                }

                for (int ires=ires1;ires<ires2;++ires) {
                    if ((ihres=ih[ires])<0) continue;

                    if (chainof[ihres]==chainof[jhres])
                        eterm=hp_pair(ihres,jhres,AtomCoordinates::s2);
                    else eterm=hp_pair(ihres,jhres,AtomCoordinates::dist2);

                    delv+=(de=eterm-Mehp[ihres*NHAA+jhres]);

                    if (de!=0) {
                        Vehp[changed[nchanges++]=ihres*NHAA+jhres]=eterm;
                        Vehp[changed[nchanges++]=jhres*NHAA+ihres]=eterm;
                    }
                }

                ++jres;
            }
        }

        return (delv*=sclfct);
    }

    void HydrophobicityPi::Accept(Update *updt)
    {
        for (;nchanges;--nchanges)
            Mehp[changed[nchanges-1]]=Vehp[changed[nchanges-1]];

        vval+=delv;
    }

    double HydrophobicityPi::hp_pair(int i, int j, double(*distf)(int,int))
    {
    	
    	//prf::cout << "hp_pair\n";
    	
	    double ans;//,cpstr;
        /* i and j are indices for NHAA sized arrays. */
	    
	    /*
	    char olc1 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[i]))->OLC() );
	    int  id1  = ((AminoAcid*) p->amino_acid(hlig[i]))->UniqueId();
	    char olc2 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[j]))->OLC() );
	    int  id2  = ((AminoAcid*) p->amino_acid(hlig[j]))->UniqueId();
	    prf::cout << "Res:"<<olc1<<id1+1 <<",Res:"<<olc2<< id2+1 << "\n";
	    prf::cout << "Pi?"<<ispi[i]<<ispi[j]<<"\n";
	    prf::cout << "Cat?"<<iscat[i]<<iscat[j]<<"\n";*/
	    bool success = false;
	    if ( (ispi[i]==true && ispi[j]==true ) || (ispi[i]==true && iscat[j]==true) || (iscat[i]==true && ispi[j]==true) ){
        //if ( validCatPi[make_pair(i,j)]){
        	//prf::cout << "Catpistrength:"<< catpistrength<< "\n";
			//ans=-strength[i*NHAA+j] * catpistrength * catpi(i,j);
        	//assert (not (iscat[i]==true && iscat[j]==true));
        	
        	/*if (iscat[i] || iscat[j])
        		
        		cpstr = 2 * catpistrength ;
        	else
        		cpstr = catpistrength;
        	
        	ans= - ( (strength[i*NHAA+j] * contact_frac(i,j,distf))  + (cpstr * catpi(i,j)) )    /2.0;*/
	    	
	        //ans = -catpistrength*catpi(i,j,success);
	    	
	    	if (mode == "1" || mode=="geometric"){
	    		ans = arpistrength*arpi(i,j,success);
	    	}else if (mode == "2" || mode=="geometricKeepPRFStrength"){
	    		ans = -strength[i*NHAA+j]*arpistrength*arpi(i,j,success);
	    	}else if (mode == "3" || mode=="avgGeometricContact"){
	    		ans = -0.5*strength[i*NHAA+j]* (arpistrength*arpi(i,j,success)  + contact_frac(i,j,distf));
	    	}
	    	//ans = -0.5*strength[i*NHAA+j]* (catpistrength*catpi(i,j,success)  + contact_frac(i,j,distf));
	        /*
	        if (ans > 0.0){
	        	ans *= cprep;
	        }
	        */
	        // when side chains are not within the range covered by the look-up table
	        if (success == false){
	        	ans=-strength[i*NHAA+j]*contact_frac(i,j,distf);
	        }
        
        	//assert(ans==0.0);
        } else {
        	ans=-strength[i*NHAA+j]*contact_frac(i,j,distf);
        }
        
        /*
        if (ans != 0.0){
                		
                		char olc1 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[i]))->OLC() );
                			    int  id1  = ((AminoAcid*) p->amino_acid(hlig[i]))->UniqueId();
                			    char olc2 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[j]))->OLC() );
                			    int  id2  = ((AminoAcid*) p->amino_acid(hlig[j]))->UniqueId();
                		
                		if (success==true) prf::cout << "IsPi! Res:"<<olc1<<id1+1 <<",Res:"<<olc2<< id2+1 << ": " << ans << " , cpstr="<<catpistrength<<", cprep="<<cprep<< "\n";
                		else prf::cout << "NotPi! Res:"<<olc1<<id1+1 <<",Res:"<<olc2<< id2+1 << ": " << ans << " , cpstr="<<catpistrength<<", cprep="<<cprep<< "\n";
                		
                }
    	*/
        
        //if (ans>0) prf::cerr <<"hp_pair: "<<i<<"\t"<<j<<" : "<<ans<<"\n";
        
        if ( ( hlig[i]==focalRes && hlig[j]==focalRes2 )
        		|| ( hlig[i]==focalRes2 && hlig[j]==focalRes ) ){
        	ans *= focalResScale;
        }
        
        if (debug){
	        if ( hlig[i] != -1 &&  hlig[j] != -1 && ans != 0.0){
	        	char olc1 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[i]))->OLC() );
    		    int  id1  = ((AminoAcid*) p->amino_acid(hlig[i]))->UniqueId();
    		    char olc2 = Groups::mapOLC2Char( ((AminoAcid*) p->amino_acid(hlig[j]))->OLC() );
    		    int  id2  = ((AminoAcid*) p->amino_acid(hlig[j]))->UniqueId();
	        		    
	        	
	        	if (debug) prf::cout << "HydrophobicityPi::hp_pair Res:"<<olc1<<id1+1 <<",Res:"<<olc2<< id2+1 << ": " << ans << "";
	        	if ( hlig[i]==focalRes || hlig[j]==focalRes){
	        		if (debug) prf::cout << " ***\n";
		        	//if ( hlig[j] != -1 && ((AminoAcid*) p->amino_acid(hlig[j]))->OLC() == F && hlig[j]==51){
		        	//	ans *= 1.0;
		        	//prf::cout << "Pi:" << hlig[i] << "," << hlig[j] << "\n";
		        
		        	//}
	        	}else if (ispi[i] && ispi[j]){
	        		if (debug) prf::cout << " pipi\n";
	        	}else{
	        		if (debug) prf::cout << "\n";
		        }
	        	
	        }
        }
        
        return ans;
    }
    
    void HydrophobicityPi::initPiInteractionMap(){
    	prf::cout << "HydrophobicityPi::initPiInteractionMap\n";
    	double tmp_r=-999.999, tmp_theta=-999.999, tmp_phi=-999.999, tmp_f;
    	double tmp_rstep=-999.999, tmp_thetastep=-999.999, tmp_phistep=-999.999;
    	
    	std::string filename, key,rkey, line;
    	
    	int expectedElements = (int) (rmax-rmin)/rstep * (thetamax-thetamin)/thetastep * (phimax-phimin)/phistep ;
    	prf::cout << expectedElements << "\n";
    	
    	for (int i=0; i< nPiFilenames; i++){
    		std::deque<std::string> lines;
    		
    		filename = pipimap_filenames[i];
    		
    		prf::cout << "Looking for file: *" << filename << "\n";
    		
    		// look for files having <filename> as a suffix
    		char* cwd = ".";
    		std::string new_fname = prf_utils::getFileContaining(cwd, filename.c_str(),true);
    		
    		bool success;
    		if (new_fname != "null"){
    			success = prf_utils::get_lines(new_fname,lines);
    		}else{
    			success = prf_utils::get_lines(filename,lines);
    		}
    			
    		
    		
    		prf::cout << lines.size()*sizeof(lines[0]) << " " << lines.max_size() / sizeof(lines[0])<< "\n";
    		//assert (lines.size() < lines.max_size());
    		
    		key = filename.substr(0,2);
    		rkey = key; // reverse key, .i.e "FY" -> "YF"
    		std::reverse(rkey.begin(), rkey.end());
    		prf::cout << "key,rkey: "<< key << ", " << rkey << "\n";
    		
    		if ( success==0) {
    			prf::cout << "Could not find:" << filename << "\n";
    			//bannedKeys.push_back(key);
    			//bannedKeys.push_back(rkey);
    			validCatPi[key] = false;
    			validCatPi[rkey] = false;
    			continue;
    		}else{
    			validCatPi[key] = true;
    			validCatPi[rkey] = true;
    		}
    		// if file does not exist: continue!
     		
    		
    		
    		
    		std::map< std::string, double>* data;
    		
    		// check if this key is already there
    		std::map<std::string, std::map<std::string, double>* >::iterator findkeyouter = AA2PiMap.find(key);
    		
    		if(findkeyouter == AA2PiMap.end()){
    			prf::cout << "added key:" << key << "\n";
    			data = new std::map< std::string, double>;
    			AA2PiMap.insert(std::make_pair(key,data));
    			pipikeys.push_back(key);
    			if (key!=rkey){
    				AA2PiMap.insert(std::make_pair(rkey, data));
    				pipikeys.push_back(rkey);
    				prf::cout << "added key:" << rkey << "\n";
    			}
    			
    		}
    		//data = AA2PiMap[key];
    		
    		
    		/*
    			std::map<std::string, std::map<std::vector<double>, double>* >::iterator findkeyouter2 = AA2PiMap.find(rkey);
    			
    			if(findkeyouter2 == AA2PiMap.end()){
    				prf::cout << "added rkey:" << rkey << "\n";
	    			
	    		}
    			//rdata = AA2PiMap[rkey];
    		*/
    		
    		// infer min, max, step for r, theta, phi directly from data file
    		rmin = -999.999;
    		thetamin = -999.999;
    		phimin = -999.999;
    		
    		rstep = -999.999;
    		thetastep = -999.999;
    		phistep = -999.999;
    		
    		int excounter = 0;
    		int counter = 0;
    		for (int j=0; j< lines.size();j++ ){
    			line = prf_utils::trim_str(lines[j]);
    			//prf::cout << "Line: "<< line << "\n";
    			std::deque<std::string> parts;
    			prf_utils::split(line,parts);
    			//prf::cout << "nParts:" << parts.size()<< "\n";
    			//double tmp_arr[];
    			if (parts.size()==4){
    				tmp_r =     abs(prf_utils::stringtodouble( parts[0] ));
    				tmp_theta = abs(prf_utils::stringtodouble( parts[1] ));
    				tmp_phi =   abs(prf_utils::stringtodouble( parts[2] ));
    				
    				
    				if (tmp_rstep != -999.999 && tmp_r != tmp_rstep)
    					rstep = tmp_r - tmp_rstep;
    				if (tmp_thetastep != -999.999 && tmp_theta != tmp_thetastep)
    					thetastep = tmp_theta - tmp_thetastep;
    				if (tmp_phistep != -999.999 && tmp_phi != tmp_phistep)
    					phistep = tmp_phi - tmp_phistep;
    				
    				tmp_rstep = tmp_r;
    				tmp_thetastep = tmp_theta;
    				tmp_phistep = tmp_phi;
    				
    				if (rmin==-999.999)
    					rmin = tmp_r;
    				if (thetamin == -999.999)
    					thetamin = tmp_theta;
    				if (phimin == -999.999)
    					phimin = tmp_phi;
    				
    				tmp_f =     prf_utils::stringtodouble( parts[3] );
    				std::ostringstream strs;
    				strs   << tmp_r << "," << tmp_theta << "," << tmp_phi;
    				std::string datakeystring = strs.str();
    				
    				
    				//prf::cout <<datakeystring<<" "<< key << ", " << tmp_r<< ", " << tmp_theta<< ", " << tmp_phi << ","<< abs(remainder(tmp_r,rstep))<<"\n";
    				//assert (tmp_r >= rmin && tmp_r <= rmax && abs(remainder(tmp_r,rstep))<0.00001);
    				//assert (tmp_theta >= thetamin && tmp_theta <= thetamax && abs(remainder(tmp_theta,thetastep))==0.0);
    				//assert (tmp_phi >= phimin && tmp_phi <= phimax && abs(remainder(tmp_phi,phistep))==0.0);
    				
    				//tmp_arr = {tmp_r,tmp_theta,tmp_phi};
    				//tmp_arr = {tmp_r,tmp_theta,tmp_phi};
    				//std::vector<double> tmp_arr;
    				//double tmp_arr[3] = {-1.0,-1.0,-1.0};
    				//tmp_arr.push_back(tmp_r);
    				//tmp_arr.push_back(tmp_theta);
    				//tmp_arr.push_back(tmp_phi);
    				std::pair<std::map<std::string,double>::iterator,bool > ret;
    				ret = AA2PiMap[key]->insert(std::pair<std::string,double>(datakeystring, tmp_f) );
    				
    				if (ret.second==false) {
    				    //std::cout << "element 'z' already existed";
    				    //std::cout << " with a key of " << ret.first->first << '\n';
    				    //std::cout << " with a value of " << ret.first->second << '\n';
    				    //prf::cout << key << ", " << tmp_r<< ", " << tmp_theta<< ", " << tmp_phi << ","<< AA2PiMap[key]->size() <<"\n";
    				    excounter++;
    				}else{
    					counter++;
    				}
    				
    				//prf::cout << key << ", " << tmp_r<< ", " << tmp_theta<< ", " << tmp_phi << ","<< AA2PiMap[key]->size() <<"\n";
     				
    			}else{
    				prf::cout << parts.size() << "\n";
    			}
    			
    		}
    		if (debug==true){
    			prf::cout << "existed:" << excounter << "\n";
    			prf::cout << "unique:" << counter << "\n";
    		}
    		
    		
    		
    		rmax = tmp_r+rstep;
    		thetamax = tmp_theta+thetastep;
    		phimax = tmp_phi+phistep;
    		
    		std::cout<<"r:"<<rmin<<","<<rmax<<","<<rstep<<"\n";
    		std::cout<<"theta:"<<thetamin<<","<<thetamax<<","<<thetastep<<"\n";
    		std::cout<<"phi:"<<phimin<<","<<phimax<<","<<phistep<<"\n";
    	}
    	
    	if (debug==true){
	    	// exhaustive lookup !!
	    	int elcounter = 0;
	    	int succounter = 0;
	    	int failcounter = 0;
	    	for (double i=rmin; i<rmax; i+=rstep){
	    		for (double j=thetamin; j<thetamax; j+=thetastep){
	    			for (double k=phimin; k<phimax; k+=phistep){
	    				elcounter++;
	    				for (int l=0; l< pipikeys.size(); l++){
	    					
	    					bool success=true;
	    					double fen = lookupPi(pipikeys[l], i, j, k, success);
	    					if (success == true) {
	    						succounter++;
	    						
	    					}else{
	    						failcounter++;
	    						//prf::cout << pipikeys[l]<< " " << i << " " << j << " " << k << " " << fen << "---";
	    					}
	    				}
					}
	    		}
	    	}
	    	
	    	prf::cout << "Done exhaustive lookup test!\n";
	    	prf::cout << "expected per file: "<< elcounter << "\n";
	    	prf::cout << "found:"<<succounter<<"\n";
	    	prf::cout << "found not:"<<failcounter<<"\n";
	    	prf::cout << "total:"<<failcounter+succounter<<"\n";
	    	//assert
	    	typedef std::map<std::string, std::map<std::string, double>* >::iterator it_type;
	    	typedef std::map<std::string, double>::iterator it_type2;
	    	for(it_type iterator = AA2PiMap.begin(); iterator != AA2PiMap.end(); iterator++) {
	    		std::map<std::string, double>* val = iterator->second;
	    		prf::cout << "key:" << iterator->first << ", len(value): " << val->size() << "\n";
	    		assert (val->size()==elcounter);
	    		//for(it_type2 iterator2 = val->begin(); iterator2 != val->end(); iterator2++) {
	    			//prf::cout <<"value: " << iterator2->second << "\n";
	    		//}
	    		
	    		// iterator->first = key
	    	    // iterator->second = value
	    	    // Repeat if you also want to iterate through the second map.
	    	}
    	}
    	//assert (false);
    }
    
    double HydrophobicityPi::getIndex(double num, double min, double factor){
    	double rem = fmod(num, factor);
    	double exact = num - rem;
    	//int index = (int) (exact-min) / factor;
    	//prf::cout << num << " " << factor << " " << rem << " " << exact << "\n";
    	
    	return exact;
    }
    
    double HydrophobicityPi::lookupPi(std::string key, double r, double theta, double phi, bool &success){
    	/*std::vector<double> values;
    	values.push_back(getIndex(r, 	rmin, 		rfactor) );
    	values.push_back(getIndex(theta, thetamin, 	thetafactor) ); 
    	values.push_back(getIndex(phi, 	phimin, 	phifactor) );*/
    	
    	//prf::cout<<"lookupPi()\n";
    	
    	double r_ind = getIndex(r, 	rmin, 		rstep);
    	double t_ind = getIndex(theta, thetamin, 	thetastep);
    	double p_ind = getIndex(phi, 	phimin, 	phistep);
    	
    	//prf::cout<<"got indices\n";
    	
    	std::ostringstream strs;
		strs  << r_ind << "," << t_ind << "," << p_ind;
		std::string datakeystring = strs.str();
    	
    	std::map<std::string, double>::iterator findkeyinner = AA2PiMap[key]->find(datakeystring);
    	
    	//prf::cout << "r:" << r << "-->" << r_ind << ", t:" << theta << "-->" << t_ind << ", p:" << phi << "-->" << p_ind << " :: "<< datakeystring<<"\n";
    	
    	//prf::cout << "found key\n";
    	if (findkeyinner == AA2PiMap[key]->end()){
    		//prf::cout << r << " " << theta << " " << phi << "\n";
    		if (debug == true)
    			prf::cout << "values NOT found! "<<r_ind<<","<<t_ind<<","<< p_ind<<" "<< datakeystring << "\n";
    		//assert (findkeyinner != AA2PiMap[key]->end());
    		success = false;
    		return 0.0;
    	}else{
    		//prf::cout << "values found! "<<values[0]<<","<<values[1]<<","<<values[2]<<"\n";
    		
    		success = true;
    		//prf::cout<<key<<", "<<r<<", "<<theta<<", "<< phi<<", "<<findkeyinner->second<<"\n";
    		return findkeyinner->second;
    	}
    	
    }

    double HydrophobicityPi::hp_contact_frac(int iaa, int jaa,
                                           double(*distf)(int,int))
    {
        /* iaa and jaa are amino acid indices relative to the
        entire system. This function is intended only for use
        in connection with hydrophobicity related observables.*/
        int i=ih[iaa],j=ih[jaa];

        if (i<0||j<0) return 0.0;

        return contact_frac(i,j,distf);
    }
    
    std::string HydrophobicityPi::key4Pair(Ligand* lig1, Ligand* lig2){
    	//prf::cout<<"key4pair:";
    	std::string tmp = "";
    	std::string s1 = "";
    	std::string s2 = "";
    	stringstream ss1,ss2;
    	ss1 << Groups::mapOLC2Char(lig1->OLC());
    	ss1 >> s1;
    	//prf::cout<<s1;
    	tmp.append(s1);
    	ss2 << Groups::mapOLC2Char(lig2->OLC());
    	ss2 >> s2;
    	tmp.append(s2);
    	//prf::cout<<s2<<"\n";
    	return tmp;
    }
    
    double HydrophobicityPi::arpi(int i, int j, bool &success){
    	/* i and j are indices for NHAA sized arrays. */
    	
    	Ligand* lig1 = p->ligand(hlig[i]);
    	Ligand* lig2 = p->ligand(hlig[j]);
    	
    	std::string lookupkey = key4Pair(lig1,lig2);
    	if (validCatPi[lookupkey] == false){
    		success = false;
    		return 0.0;
    	}
    	
    	Vector3 c1,cat1a0,cat1a1;
    	if (ispi[i])
    		c1 = centroid(lig1);
    	else if (iscat[i]){
    		cat1a0 = AtomCoordinates::vec(lig1->labeled_atom(catatms[ctyp(lig1->OLC())][0])->UniqueId());
    		cat1a1 = AtomCoordinates::vec(lig1->labeled_atom(catatms[ctyp(lig1->OLC())][1])->UniqueId());
    		c1 = cat1a0;
    	}
    	
    	Vector3 c2,cat2a0,cat2a1;;
    	if (ispi[j])
    		c2 = centroid(lig2);
    	else if (iscat[j]){
    		cat2a0 = AtomCoordinates::vec(lig2->labeled_atom(catatms[ctyp(lig2->OLC())][0])->UniqueId());
    		cat2a1 = AtomCoordinates::vec(lig2->labeled_atom(catatms[ctyp(lig2->OLC())][1])->UniqueId());
    		c2 = cat2a0;
    	}
    	
    	Vector3 cc = c2-c1;
    	double r = cc.mag();
    	
    	if (r < rmin || r > rmax){
    		return 0.0;
    	}
    	
    	Vector3 n1;
    	if (ispi[i])
    		n1 = normal(lig1);
    	else if (iscat[i]){
    		n1 = cat1a1 - cat1a0;
    	}
    	
    	Vector3 n2;
    	if (ispi[j])
    		n2 = normal(lig2);
    	else if (iscat[j]){
    		n2 =  cat2a1 - cat2a0;
    	}

    	double theta = n1.angle(n2) ;
    	theta = min(theta, M_PI-theta)* (180.0/M_PI);
    	
    	if ( (theta < thetamin) || (theta > thetamax)){
    		return 0.0;
    	}
    	
    	double phi = n1.angle(cc) ;
    	phi = min(phi, M_PI-phi) * (180.0/M_PI);
    	
    	if ( (phi < phimin) || (phi > phimax) ){
    		return 0.0;
    	}
    	
    	success = true;
    	double value = lookupPi(lookupkey, r, theta, phi, success);
    	
    	if (debug==true)
    		if (value!=0.0)
    			prf::cout << "r, theta, phi, F: " << r << " " << theta << " " << phi << " " << value <<"\n";
    	assert (success);
    	return value;
    }
    
    double HydrophobicityPi::contact_frac(int i, int j,
                                        double(*distf)(int,int))
    {
        /* i and j are indices for NHAA sized arrays. */
        int ni,nj,k,l,m,n;
        double tmp,sum=0;
        ni=natom[i];
        nj=natom[j];

        for (n=0;n<ni+nj;n++) r2min[n]=1e10;

        for (n=0;n<ni;n++) {
            k=atom[i*6+n];

            for (m=0;m<nj;m++) {
                l=atom[j*6+m];
                tmp=distf(k,l); // distance between atoms k,l

                if (tmp>max_at_dist2 && tmp<max_cmp_dst2) return 0;

                if (tmp<r2min[n]) r2min[n]=tmp;

                if (tmp<r2min[ni+m]) r2min[ni+m]=tmp;
            }
        }

        for (n=0;n<ni+nj;n++) {
            
        	if (r2min[n]>b2) continue;
        	//prf::cout << r2min[n] << "\n";
            if (r2min[n]<a2) sum++;
            else sum+= (b2-r2min[n]) /d2;
        }

        double gammaij=1;

        if ((soften[i] == soften[j]) and (soften[i]!=0)) gammaij=0.75;
        
        //prf::cout << "contact_frac=" << std::min(1.0,sum/(gammaij*(ni+nj))) << ", " << "\n";
        return std::min(1.0,sum/(gammaij*(ni+nj)));
    }
}

double HydrophobicityPi::dhp_pair(int i, int j, bool perdhint,
                                std::valarray<double> &gx)
{
    double frac=0,pstr=0;
    int ni,nj,k,l,m,n;
    double tmp,sum=0;
    ni=natom[i];
    nj=natom[j];
    pstr=-strength[i*NHAA+j]*sclfct;
    for (n=0;n<ni+nj;++n) r2min[n]=1e10;

    for (n=0;n<ni;++n) {
        k=atom[i*6+n];

        for (m=0;m<nj;++m) {
            l=atom[j*6+m];
            tmp=perdhint?AtomCoordinates::dist2(k,l):AtomCoordinates::s2(k,l);

            if (tmp>max_at_dist2 && tmp<max_cmp_dst2) return 0;

            if (tmp<r2min[n]) {
                r2min[n]=tmp;
                r2minpartner[n]=m;
            }

            if (tmp<r2min[ni+m]) {
                r2min[ni+m]=tmp;
                r2minpartner[ni+m]=n;
            }
        }
    }

    for (n=0;n<ni+nj;++n) {
        if (r2min[n]>b2) continue;

        if (r2min[n]<a2) sum++;
        else sum+= (b2-r2min[n]) /d2;
    }

    double gammaij=1;

    if ((soften[i] == soften[j]) and(soften[i]!=0)) gammaij=0.75;

    frac=(gammaij*(ni+nj));
    if (sum<=0) return 0;
    else if (sum<frac) {
        pstr/=(frac);
        for (n=0;n<ni+nj;++n) {
            if (r2min[n]>b2 or r2min[n]<a2) continue;
            if (n<ni) {
                m=atom[i*6+n];
                l=atom[6*j+r2minpartner[n]];
            } else {
                m=atom[6*j+n-ni];
                l=atom[6*i+r2minpartner[n]];
            }
            Vector3 dirn;
            if (perdhint) dirn=AtomCoordinates::sep(m,l);
            else dirn=AtomCoordinates::diff(m,l);
            dirn*=(2*pstr*slpe);
            gx[3*m]+=dirn.x();
            gx[3*m+1]+=dirn.y();
            gx[3*m+2]+=dirn.z();
            gx[3*l]-=dirn.x();
            gx[3*l+1]-=dirn.y();
            gx[3*l+2]-=dirn.z();
        }
        return pstr*sum;
    } else return pstr;
}
