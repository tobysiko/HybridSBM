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

#ifndef UPDATESHANDLER_HH
#define UPDATESHANDLER_HH
#include "../Aux/HandlerBase.hh"
#include "../Aux/RandomNumberBase.hh"
#include "Rot.hh"
#include "BGS.hh"
#include "Pivot.hh"
#include "Rotation.hh"
#include "Translation.hh"
#include <vector>

//! A manager for MC updates
/**
  * UpdatesHandler manages a lot of things relating to conformational
  * updates for Monte Carlo simulations. It's the interface for the
  * conformational updates module in ProFASi. Based on a given Population,
  * it can automatically set up a reasonable choice of conformational
  * updates. It can calculate good relative frequency of performing those
  * updates. The UpdateHandler is also used as a helper class in MC to
  * bring about the conformational changes in each MC step.
  *
  * The behaviour of the update manager can be controlled through commands
  * passed to it as InstructionStrings. For instance, it can be told to
  * use an update although the automatic set up does not select that
  * update using its more generic criteria.
  *
  */
class UpdatesHandler : public HandlerBase
{
public:
    //! Default constructor
    UpdatesHandler();
    ~UpdatesHandler();
    //! Parse updates related instructions
    int parseCommand(InstructionString s);
    //! Smart choice of updates
    /**
      * For single chain systems, nothing is normally gained by performing
      * rigid body updates. So, this function checks if there is indeed
      * only one chain, and switches off those updates... that is, if the
      * user did not explicitly state that such an update is desired! For
      * instance, the single protein might be subject to an external force
      * field which is space dependent, in which case, even rigid body
      * updates are dynamically interesting. In such a situation, the user
      * should make an explicit call to useUpdate("Rotation"), for instance,
      * which would tell the autoSelectUpdates() function to leave the
      * rigid body rotations alone, even if there is only one chain.
      */
    void autoSelect(Population *p);
    //! Number of updates selected for the given Population
    inline size_t num_updates() const {return used.size();}
    //! Get a pointer to the i'th update in use
    Update * used_update(size_t i) { return (i<used.size())?used[i]:NULL; }
    //! Get a pointer to an update by name, if it is in use
    Update * used_update(std::string upnm);
    //! Whether an update of a given name is in the list for the given Population
    bool have_update(std::string upnm);
    //! Explicitly state that an update should be used
    int use_update(Update * updt);
    //! Explicitly mark one of the known updates by name, to be used.
    int use_update(std::string updtname);
    //! Request that one model update be skipped
    int skip_update(std::string updtname);
    bool is_update_name(std::string upnm);
    //! Propagate temperature index range to updates
    void set_n_temps(size_t nt);
    //! Assign probabilities from file or use auto assignment
    /**
      * This function first tries to read in update probabilities
      * from a file in the format described in \ref update_probs.
      * If that fails, it uses population information to set up some
      * reasonable values.
      */
    void assign_probs();
    inline double operator()(size_t iup,size_t jtmp) const {
        return probs[iup][jtmp];
    }
    void set_population(Population *popl);
    void setBeta(double bt);
    void print_updateprobs(Output &op);
    void print_setup();
    void output_statistics(std::string flnm);
    void reset_statistics();
    inline void RandomNumberGenerator(RandomNumberBase *ran) {
        rng=ran;
        for (size_t i=0;i<used.size();++i) used[i]->set_RandomNumberGenerator(rng);
    }
    //! Perform any necessary initialisation operations for the updates
    void init();
    Update *perform_update(int itmp=0);
    void accept_update();
    void reject_update();
    Update * new_update(std::string upnm);
private:
    //! Set up probabilities for different updates from dof data
    /**
     * The function analyzes the system to find the number of dof of
     * different kinds: backbone, side chain, rigid body.. It then assigns
     * proportionate values for the probabilities for updates involving
     * backbone, side chain and rigid body degrees of freedom.
     */
    void auto_assign_probs();
    //! Read probabilities from a file
    /**
    * The syntax for specifying probabilities is given in \ref update_probs
    */
    int ReadFile(std::string filename);
    std::vector<Update *> used;
    void normalize();
    Matrix<double> probs;
    Matrix<size_t> num_calls, num_acc_calls;

    std::string updtprobfile;
    size_t ntmp,uplog;
    int active_update,tindex;
    RandomNumberBase *rng;
    Population *p;
};

#endif // UPDATESHANDLER_HH
