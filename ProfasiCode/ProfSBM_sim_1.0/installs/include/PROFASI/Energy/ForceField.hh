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

#ifndef FORCEFIELD_HH
#define FORCEFIELD_HH
#include <string>
#include <vector>
#include "Energy.hh"
#include "../Updates/Update.hh"

namespace prf {
    class ForceField
    {
    public:
        ForceField();
        ForceField(std::string ffname);
        ~ForceField();
        inline void set_name(std::string nm) {myname=nm;}
        inline std::string Name() const { return myname; }
        inline size_t n_terms() const { return e.size(); }
        inline Energy * term(int i) { return e[i]; }
        inline void add_term(Energy *et) {
            e.push_back(et); incharge.push_back(et);
        }
        inline void add_external_term(Energy *et) {e.push_back(et);}
        int delete_term(std::string nm);
        inline void connect(prf::Population *popl) {
            for (size_t i=0;i<e.size();++i) if (e[i]!=NULL) e[i]->Connect(popl);
        }
        void init();
        inline void set_Tmix(double t){
        	Tmix=t;
        	doMix = (Tmix != 0);
        	}
        inline double get_Tmix(){return Tmix;}
        inline void set_scaleSBM(double n){scaleSBM=n;}
        inline double get_scaleSBM(){return scaleSBM;}
        inline void set_arpi(double n){arpi=n;}
        inline double get_arpi(){return arpi;}
                
        Energy * term_called(std::string ename);
        std::string summary();
        void print_contributions(prf::Output &op);
        void refresh();
        inline double value() { return etot; }
        inline double deltaE() { return dele; }
        double evaluate();
        double evaluateLambda(double lambda);
        double deltaE(Update *updt);
        double deltaE(Update *updt, double maxde);
        void accept(Update *updt);
        void reject(Update *updt);

        double reset_total();
        double reset_total_silently();
        inline void setMixable(std::vector<Energy *> m){mixable = m;}
        inline std::vector<Energy *> getMixable(){return mixable;}
        inline void setDebug(bool b){debugFF = b;}
        inline bool hasNatCon(){return hasNative;}
        inline double getNatEn(){return SBM_E;}
        double SBM_Energy();
        double SBM_EnergyLambda(double lambda);
        double SBM_Evaluate();
        double SBM_DeltaE(Update *updt);
        double SBM_Minimum();
        inline void set_SbmToMix(std::string s){
        	s2m.resize(s.length());
        	for (size_t i =0; i < s.length(); i++){
        		if (s.at(i)=='1')
        			s2m[i] = true;
        		else
        			s2m[i] = false;
        		std::cout<<"s2m["<<i<<"]="<<s2m[i]<<"\n";
        	}
        }
        double Mix(double d);
        double scaledMix(double d, double eminsum);
        inline void set_mixtype(int m){mixtype = m;}
    private:
    	bool hasNative, debugFF, doMix;
    	std::vector<bool> s2m;
        std::vector<Energy *> e,incharge;
        std::string myname;
        double dele,etot,SBM_E,SBM_E_prop, Mlow,Mhigh;
        std::vector<double> SBM_energies_tmp,sigconsts,minE;
        std::vector<Energy *> mixable;
        double Tmix,scaleSBM,arpi,epsilon,nRestraints,sigconst_mix;
        int mixtype;
    };

}

class Etot : public Observable
{
public:
    Etot();
    ~Etot();
    inline void connect(ForceField *ff) { myff=ff;}
    double evaluate();
    double delta(Update *u);
    void refresh();
    void rangeEstimate(double &x1,double &x2);
private:
    ForceField *myff;
};

class Ego : public Observable
    {
    public:
    	Ego();
        ~Ego();
        inline void connect(ForceField *ff) { myff=ff;}
        double evaluate();
        double delta(Update *u);
        void refresh();
        void rangeEstimate(double &x1,double &x2);
    private:
        ForceField *myff;
    };

#endif // FORCEFIELD_HH
