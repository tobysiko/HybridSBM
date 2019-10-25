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

#ifndef OBS_ENERGY_HH
#define OBS_ENERGY_HH
#include "../Energy.hh"

namespace prf
{
    //! This class provides an energy term based on an arbritrary Observable
    /**
    * It is sometimes desirable to be able to use a measurement performed on a
    * structure as a "constraint", or an energy term driving the dynamics towards
    * a larger or smaller values of that Observable. This class helps with one
    * approach towards such constrained simulations. Any Observable, such as
    * RMSD, native contacts, radius of gyration etc. can be used to guide the
    * simulations, through this pseudo-energy term.
    *
    * \sa obs_en_restraints
    */

    class ObsEnergy : public Energy
    {
    public:
        //! Default constructor
        ObsEnergy();
        ~ObsEnergy();
        //! Set scale
        void setScaleFactor(double x);
        //! Get scale factor
        inline double getScaleFactor() const {return sclf;}
        //! Connect to one Observable
        void connectObs(Observable* o);
        double evaluate();
        double gradientXYZ(std::valarray<double> &g);
        double gradientDOF(std::valarray<double> &ans,std::vector<int> &indxs);
        //! Name of the observable associated with the energy term
        std::string obsName() const ;
        //! Set name of the observable this term should be associated with
        void setObsName(std::string s);
        //! Estimate a range
        void rangeEstimate(double &x1, double &x2);
        //! Configure the ObsEnergy term
        /**
          The term can be configured by a string consisting of a comma separated
          pair of assignments like this: "obs=rg,scale=50".
          */
        void set_pars(std::string pars);
    private:
        double sclf;
        Observable* obs;
        std::string obs_name;
    };
}
/**
  \page obs_en_restraints Using an arbitrary measurement as a restraint

  ProFASi measurements module provides a large number of measurements such
  as RMSD, radius of gyration, secondary structure etc. If desired, any of
  those measurements can be used as a pseudo-energy term to restrain the
  simulation. While this does not make much sense for a folding simulation
  of one protein, in many studies  some form of restraints are useful. For
  instance, if we want to study the binding of one small peptide with a
  large protein, large scale conformational changes of the big protein are
  not the prime objects of the study. A free MC simulation will normally
  unfold the big protein and try to explore its entire conformation space --
  that's what an MC run is supposed to do, i.e., explore conformation space
  without being trapped in any minimum for ever. One way to keep a protein
  near a given structure is to use the RMSD of the current state with the
  reference state as an energy term. Larger the RMSD, larger the penalty.

  Instead of anticipating all possible measurements which someone might want
  as an energy term in a study, we provide a generic mechanism to use
  just about any measurement done during a run as an energy term. The syntax
  looks like this in the settings file:

  \verbatim
  force_field FF08Rg=FF08+Extras:ObsEnergy(obs=rg,scale=50)
  \endverbatim
  This sets up a force field that is the sum of the default force field and
  a term that is 50 times the value of an observable called "rg". Presumably
  the settings file also has a line creating a measurement called "rg":
  \verbatim
  new_obs Rg rg of_chain 0
  \endverbatim
  The name "rg" used is the alias given to the measurement by the user, not the
  name of the class. Note that the above scale factor of 50 on the radius of
  gyration is just a random example. With this value, almost any protein with
  30 residue or more is likely to turn into a nice ball with hardly any
  secondary structure.

 \warning Energy terms themselves are measurements in ProFASi. So, it is possible
to generate nonsense by "connecting" the new ObsEnergy object to itself!

 */
#endif
