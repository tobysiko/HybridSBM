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

#ifndef FFHANDLER_HH
#define FFHANDLER_HH
#include "../Elements/Population.hh"
#include "../Aux/HandlerBase.hh"
#include "ForceField.hh"
#include "Extras/ObsEnergy.hh"
#include <string>
#include <list>
#include <sstream>

class FFHandler : public HandlerBase
{
public:
    FFHandler();
    ~FFHandler();
    int parseCommand(InstructionString s);
    int init_ff();
    int set_ff(std::string ffn);
    inline ForceField *interaction_potential() { return ff; }
    inline void set_population(prf::Population *popl) {
        p=popl;
        if (ff!=NULL) ff->connect(p);
    }
    //! Add a new energy term
    void useEnergy(Energy * en);
    //! Request that one model energy term be skipped.
    /**
         * This may be useful if you have an alternative version of the model
         * energy term that you want to try out...
         */
    void skipEnergy(std::string enname);
    void export_options_to(prf_utils::ProgArgs &pars);
    size_t n_obs_dependent_terms() {return oe.size();}
    ObsEnergy * obs_dependent_term(size_t i) {return oe[i];}
    Energy * new_energy_term(std::string ename);
    
private:
    ForceField * set_up_force_field(std::string ffname);
    bool get_closure(std::list<std::string> &lst,
                     std::list<std::string>::iterator st,
                     std::list<std::string>::iterator &nd);
    void tokenize(std::string req, std::list<std::string> &tokens);
    ForceField *ff;
    std::deque<ObsEnergy *> oe;
    Population *p;
    std::string ffname, posmap;
    double Tmix;
    bool alwaysMix;
    double scaleSBM;
    double arpi,focalResScale,cprep;
    double coop;
    std::string ppmode;
    std::string s2m;
    int focalRes,focalRes2,mixtype;
    std::vector<std::string> skippedTerms;
    std::vector<Energy *> mixable;
    
    inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }


    inline std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, elems);
        return elems;
    }

};

#endif // FFHANDLER_HH
