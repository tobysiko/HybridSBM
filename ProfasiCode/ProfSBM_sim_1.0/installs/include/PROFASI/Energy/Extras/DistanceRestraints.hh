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

#ifndef DISTANCERESTRAINTS_HH
#define DISTANCERESTRAINTS_HH
#include "../Energy.hh"
#include "RestraintFunction.hh"

namespace prf {
    //! A single distance restraint
    /**
      A distance restraint consists of a pair of atoms and a function
      of the distance between them.
      */
    class DistanceRestraint {
    public:
        DistanceRestraint();
        ~DistanceRestraint();
        void delete_restraint_function();
        DistanceRestraint(const DistanceRestraint &dr);
        DistanceRestraint &operator=(const DistanceRestraint &dr);
        int set_pars(prf_xml::XML_Node *pars, Population *p);
        double evaluate();
        double estimate_upper_bound();
    private:
        int get_aid(std::string atm, Population *p);
        int atom1,atom2;
        RestraintFunction *f;
    };

    //! Distance restraints
    /**
      This is an extra energy term which can be used to encourage the
      formation of structures with particular pairs of atoms at given
      distances. The choice of atoms as well as the functional form of
      the distance dependence is fairly general. How to use this class
      to impose distance restraints in Monte Carlo simulations is described
      in \ref ds_restraints .
      \sa \ref ds_restraints
      */
    class DistanceRestraints : public Energy
    {
    public:
        //! Default constructor
        DistanceRestraints();
        ~DistanceRestraints();
        //! Set parameters from an XML file
        void set_pars(std::string xmlfilename);
        void init();
        double evaluate();
        void rangeEstimate(double &x1,double &x2);
    private:
        std::vector<DistanceRestraint> c;
        std::string filename;
    };
}
/**
  \page ds_restraints Simulations with distance restraints
Using distance restraints in a simulation follows exactly the same syntax
as does the use of dihedral restraints (See \ref dh_restraints ). The following
one line in the settings file will set up a new energy term according to the
restraints specified in the XML configuration file \tt dist_restr.xml .
\verbatim
force_field FF08DR=FF08+Extras:DistanceRestraints(dist_restr.xml)
\endverbatim
The procedure for specifying the restraints in the command line is the same as
for the dihedral restraints.

\section specs Configuring distance restraints using an XML file
The distance restraints energy term must be configured so that it knows which
distances to use and what to do with them. The XML configuration file should
have the following format:

\li The top level node for this XML configuration file is called
      \tt distance_restraints .
\li Within the scope of the top level node, there should be a series of nodes
  of name \tt restraint .
\li Each \tt restraint node must have a field called \tt atom1 , a field
  called \tt atom2 and a field called \tt parameters
\li The \tt restraint node can optionally have an attribute \tt type,
specifying the kind of function to be used for the restraint. The default
function for distance restraints is a quadratic form.
\li The contents of the \tt parameters field depend on the type of function
chosen with the type attribute. The parameter names for the default function
are \tt mean and weight.

\section example Example
The following XML file configures 5 distance restraints:
\verbatim
<distance_restraints>

<formatted_data>

<format name="restraint" type="$3">
  <atom1>$1</atom1>
  <atom2>$2</atom2>
  <parameters>
    <mean>$4</mean>
    <weight>$5</weight>
  </parameters>
</format>

<data>
0/0/GLY/_CA_ 0/15/GLU/_CA_ quadratic 5.5 3.0
0/1/GLU/_CA_ 0/14/THR/_CA_ quadratic 5.5 3.0
0/2/TRP/_CA_ 0/13/VAL/_CA_ quadratic 5.5 3.0
0/3/THR/_CA_ 0/12/THR/_CA_ quadratic 5.5 3.0
0/4/TYR/_CA_ 0/11/PHE/_CA_ quadratic 5.5 3.0
</data>

</formatted_data>

</distance_restraints>
\endverbatim

In connection with dihedral restraints, we introduced the template interpretation
capacity of the ProFASi XML module. The above is another example. Here the
tabular data to be re-cast as a series of XML nodes appears inline inside a
special tag \tt data. The \tt format tag is used to interpret the contents
of the data tag. Here we also introduce how to specify the functional form
by giving a type attribute to the auto-generated restraint nodes.

The atoms are specified through their chain, residue number, residue type, and
atom label. The labels are 4 characters wide. Spaces in the label should be
replaced by underscores. It is important to have the labels correct. It is
" CA " (and therefore "_CA_" above), and not "CA__" for C-alpha.

Several alternative forms are available for the functional form. In place of
quadratic, one could write "power_law", "flattened_power_law" or "gaussian" to
use the functions represented by prf::PowerLawRestraint, prf::FlattenedPL and
prf::GaussianRestraint respectively.

\sa \ref dh_restraints
  */
#endif // DISTANCERESTRAINTS_HH
