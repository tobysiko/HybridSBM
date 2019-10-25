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

#ifndef Update_HH
#define Update_HH
#include "../Aux/Named.hh"
#include "../Aux/IndexSelector.hh"
#include "../Aux/prf_xml.hh"
#include "../Aux/InstructionString.hh"
#include "../Elements/Population.hh"

/**
* \defgroup profasi_updates Conformational Updates
* @brief Conformational updates for proteins and related utilities
*
* These classes represent different ways to propose a new conformation
* of a protein chain during a Monte Carlo simulation. An update could be
* a single angle update, representing a certain kind of change in one of
* the degrees of freedom, or a concerted update of a series of angles.
* An update class should conform to the interface defined by the base
* class called Update. Any class that inherits from Update and provides
* sensible actions for the different virtual member functions of Update
* (by overriding the virtual functions) can be used in a Monte Carlo
* simulation with ProFASi. In fact, it is quite likely that if an update
* class provides the information that the virtual functions in the base
* class ask for, the energy classes of ProFASi will automatically know
* how to evaluate increments in energy due to the update in an optimized
* manner. This module also contains other entities related to conformational
* updates, such as probability calculators for such updates.
*/

namespace prf
{
    //! A structure to hold the most basic information about one change
    typedef struct {
        DOF_Info info;
        double before,after;
    } dof_change_type;

    //! The base class for all conformational updates
    /**
     * In this class we define the interface that every conformational update in
     * PROFASI must have. The individual updates, like bgs for instance, would
     * over-ride the virtual functions, and define exactly what happens when that
     * update is invoked. But having a common interface facillitates integration
     * with the rest of the package. For instance, the MC class makes a vector
     * of Update's and calls the functions perform(), accept() or reject() on
     * the Update object.
     *
     * In other words, Update is an abstract concept... something happens to the
     * system conformation. That is all that is required to define the behaviour
     * of the Markov Chain generators. So, all such classes which only require
     * an abstract interface, make use of this base class, and invoke Update
     * objects through pointers. A real application on the other hand, will not
     * define abstract Updates but rather objects of the derived classes. But a
     * pointer to an Update can also point to an object of a derived class. This
     * way, the implementation of the MC class or its derivatives can be kept
     * separate from the details and finer points of individual conformational
     * updates.
     * \ingroup profasi_updates
     */

    class Update : public Named
    {
    public:
        Update();
        virtual ~Update();
        //! An update needs a random number generator
        void set_RandomNumberGenerator(RandomNumberBase * rn);
        //! Sweep through the DOFs sequencially rather than randomly
        inline void set_sequencial_sweep(bool sw=true) { issequencial=sw; }
        //! An update must be connected with a population
        void connect(Population * pl);
        //! Initialize update for use
        virtual void init();

        //! Make a move in the conformation space
        inline int perform(int itm) {itmp=itm; return perform();}
        virtual int perform();

        //! Finalizes a proposed update
        virtual int accept();

        //! Reverts the proposed update
        virtual int revert();

        //! Print info about any departure from default behaviour
        virtual void print_setup(std::string &st);

        //@{
        /**
          @name Modulating the behaviour of updates
          */
        //! The number of DOF in the system which may be touched by the update
        inline size_t n_relevant_dof() const { return site_weight.size(); }
        //! Weight of a particular DOF
        /**
          The input argument is the unique id of a DOF. The output is 0 if the
          DOF will never be touched: like if the DOF is a sidechain angle and
          the update is a purely backbone update. Otherwise, it is the relative
          weight of the DOF compared to other DOFs in the domain of the update.
          Appart from a normalising factor, this corresponds to the number of
          times a certain DOF is updated in a sweep.
          */
        double get_weight_of_dof(int idof);
        //! Set weight of a DOF
        /**
          If the DOF is in the domain of an update, this changes the probability
          that it is selected for a move. All DOFs in the domain are given a
          weight 1 by default, which means that they will be updated with
          equal probability. If a DOF is given a higher weight, it will be
          favoured. Normalisation is taken care off internally. Like get_weight_of_dof,
          the input argument is the global unique id of a DOF.
          */
        void set_weight_of_dof(int idof, double vl);
        //! Configure update using options provided in an XML node
        virtual void configure(prf_xml::XML_Node *nd);
        //! get site index of a given DOF
        /**
          Search and return the site index of a certain DOF. If the DOF is
          not in the domain of the update, the size of the dof_at_site
          array is returned.
          */
        size_t get_site_of_dof(DOF_Info &dof);

        //@}

        //@{
        /**
          @name Update categories
          Functions returning yes/no answers about if the updates is a ...
          */
        inline bool sidechain_update() const {return isscu;}

        inline bool rigid_chain_update() const {return isrgu;}

        inline bool backbone_update() const {return isbbu;}

        inline bool multichain_update() const {return ismcu;}

        inline bool local_update() const {return islocal;}
        //@}

        //@{
        /**
          @name So, what changed ?
          Functions returning information about what degrees of freedom
          changed because of an update. */
        inline unsigned num_changes() { return nchanges; }

        inline dof_change_type & change(unsigned i) { return thechange[i%nchanges]; }

        //@}

        //@{
        /**
          @name Ranges of changed things
          */
        //! The first and last atoms affected by the update
        inline void current_atoms(int &istr, int &iend) {
            istr = st_atom;
            iend = nd_atom;
        }

        //! First affected atom
        inline int begin_atom() { return st_atom; }

        //! One past the last atom, like the end() in standard containers
        inline int end_atom() { return nd_atom; }

        //! Start of the flexible part of the update.

        /**
         * If the update moves a section of the system in which lots of atoms
         * move non-rigidly wrt eachother, we define a "flexible_part".
         */
        inline int begin_flexible_part() { return st_fl; }

        //! End of the flexible part, once again, one past the last ...
        inline int end_flexible_part() { return nd_fl; }
        //! Additional weight for Metropolis-Hastings updates
        /**
         PROFASI provides a method to use a corrective weight for an
        update on top of what would be used in a pure Metropolis step.
        Such a weight is required if the probabilities for the size of
        the update somehow depend on the starting conformation, so that
        the forward and backward processes would have would have different
        distributions. In order to ensure detail balance, one needs to
        introduce an additional weight that compensates for the above.
        In PROFASI, only BGS has that kind of a weight. But the provision
        is made more generally to smoothen the interface and allow for
        more updates like the BGS in the future. The corrective weight is
        called IntrinsicWeight.
         */
        virtual double intrinsic_weight() const;

        //! Makes list of ranges of changed residues
        /**
        * Each range returned represents a set of residues which have
        * moved in a rigid manner relative to each other.
        */
        inline int n_residue_rigid_ranges() const {return nranges;}

        //! i'th rigid range returned as ligand indices in r1 and r2
        /**
        An update can divide the system into many segments which move
        rigidly with respect to the rest of the system. Since the degrees
        of freedom in PROFASI are only the torsional angles, it is most
        often the case that there are indeed only a few subsystems that
        move rigidly with respect to the rest. Most easy to visualize
        is a rigid body translation of a single chain relative to the
        rest. But also when a torsional angle is changed, two parts of
        the system move rigidly with respect to each other about the axis
        of that torsional angle.

        Here, the ranges returned are not the atomic ranges, but rather
        at the level of residues, or ligands, relative to the entire system.

        The following convention is
        used here, and is important in the design of new updates: If a
        range contains more than one ligand, atoms in one ligand have a
        fixed distance with respect to all atoms in the other ligands
        covered in the range. If the range contains only one ligand, no
        such restriction is assumed. That is, atoms inside that ligand
        can change distance between each other. Such single ligand
        ranges must be used for residues in which one or more degrees
        of freedom have changed. Longer ranges can only be used for
        residues in which no internal degrees of freedom have changed.

        The ranges are inclusive at both the left and the right, meaning,
        a range of (2,5) would include residues 2,3,4,5 as a rigidly moving
        block.

        It is further assumed that the ranges are arranged in increasing
        ligand indices, and are mutually non-overlapping, and are non-degenerate,
        meaning they always contain at-least one residue.
        */
        void residue_rigid_range(int i, int &r1,int &r2);

        //! Pointer to the vector of rigid ranges.
        inline std::vector<std::pair<int,int> > * residue_rigid_ranges() {
            return &ranges;
        }

        inline int n_atom_rigid_ranges() const {return natomranges;}

        void atom_rigid_range(int i, int &r1, int &r2);
        inline std::vector<std::pair<int,int> > *atom_rigid_ranges() {
            return &atom_ranges;
        }
        //@}

        virtual void set(InstructionString cmd);
        virtual void setScale(double gscl);

    protected:
        // Make a list of degrees of freedom controlled by an update
        virtual void build_dof_list();
        Population *popl;
        RandomNumberBase *rnd;
        IndexSelector selector;
        std::vector<dof_change_type> thechange;
        // dof_at_site is a map of dofs under the jurisdiction of an update.
        // For instance, for pivot, it will be a vector containing backbone
        // angles which can be changed, in the entire population. Each element
        // of the array contains redundant information to uniquely identify
        // a degree of freedom in many different ways.
        std::vector<DOF_Info> dof_at_site;
        // site weight is used to set up the index selector
        std::vector<double> site_weight;
        bool isscu, isrgu, isbbu, ismcu, islocal,issequencial;
        std::vector<std::pair<int,int> > ranges, atom_ranges;
        int st_atom, nd_atom, nranges, natomranges;
        int state, st_fl, nd_fl, site, itmp;
        unsigned nchanges;
    };
}

/**
  \page config_update XML configuration file for updates
  ProFASi simulation programs offer a way to configure the behaviour of
  conformational updates to a large extent using an XML configuration file. The
  configuration file is passed to the program with a command in the settings file
  or the command line as follows:
  \verbatim
  config_updates updates.xml
  \endverbatim
  \section upxml Structure of the update configuration XML file
  The top level node should be called \tt config_updates, and it should
  contain child nodes containing sections for any update you want to configure.
  The tags corresponding to the update specific nodes should be the name of the
  individual update, or "all_updates". Here is an example where we increase
  the probability of choosing backbone angles in residue indices 3,4,5 and 6 by
  a factor 2.
  \verbatim
  <config_updates>
    <Pivot>
      <dof id="0:3:b:0"> <weight> 2 </weight> </dof>
      <dof id="0:3:b:1"> <weight> 2 </weight> </dof>
      <dof id="0:4:b:0"> <weight> 2 </weight> </dof>
      <dof id="0:4:b:1"> <weight> 2 </weight> </dof>
      <dof id="0:5:b:0"> <weight> 2 </weight> </dof>
      <dof id="0:5:b:1"> <weight> 2 </weight> </dof>
      <dof id="0:6:b:0"> <weight> 2 </weight> </dof>
      <dof id="0:6:b:1"> <weight> 2 </weight> </dof>
    </Pivot>
  </config_updates>
  \endverbatim

  The degrees of freedom for which the behaviour should be changed is
  identified using the \ref dof_strings . The DOF specific nodes above can
  also be obtained elegantly from a tabular data file using ProFASi's
  XML template handling capacity (See \ref xml_templates ).
  */
#endif
