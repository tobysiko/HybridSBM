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

#ifndef POPULATIONHANDLER_HH
#define POPULATIONHANDLER_HH
#include "Population.hh"
#include "../Aux/HandlerBase.hh"
#include <list>

//! Helper class for configuring the protein population
/**
  The population handler owns a population object, and provides
  a broader interface to it. The idea is to reduce the Population
  class to a simple container of proteins, and provide well defined
  catetories of functionality in small helper classes like this one.
  At the present, PopulationHandler is only responsible for
  initializing the population: either from a ProFASi XML structure
  file, or (through a set of settings file commands) from sequence
  alone, or from a PDB file.
  */
class PopulationHandler : public HandlerBase
{
public:
    //! Default constructor
    PopulationHandler();
    ~PopulationHandler();
    //! Retrieve pointer to population
    inline Population *population() { return &p; }
    //! Set up a population from an XML population description file
    int read_xml_pop(std::string xmlfile);
    //! Set up a population from a PDB file and selection
    /**
      * This function will import a selection from a given PDB file as the
      * population. This means, any pre-existing chains in the population
      * will be lost. The selections are given in the usual method
      * \ref prf_sel_fils .
      */
    int read_pdb_pop(std::string pdbfilename, std::string sel);
    //! Parse commands passed as a deque of InstructionString objects
    int parseCommand(InstructionString s);
    int init_pop();
    int init_coords();
    int reconstruct();
    //! Set size of the periodic box
    /**
      * PROFASI simulations are always done in a periodic box. Here we
      * set the box size. If the user does not provide a box size, a sensible
      * default is calculated by analyzing the population. The parameter
      * usbit is used to distinguish a box size set up promted by a user
      * choice (settings file etc) from an automatic calculation for the
      * case when such a choice is not given. Values specified by
      * the user have priority.
      */
    void boxSize(double xx, bool usbit=true);
    void setupPeriodicBox();
private:
    prf::Population p;
    std::string initstate;
    double boxl;
    bool usrspc_box_size;
};

#endif // POPULATIONHANDLER_HH
