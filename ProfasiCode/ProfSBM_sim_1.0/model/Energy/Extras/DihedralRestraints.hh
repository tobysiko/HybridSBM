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

#ifndef DIHEDRALRESTRAINTS_HH
#define DIHEDRALRESTRAINTS_HH
#include "../Energy.hh"
#include "RestraintFunction.hh"
#include <map>

namespace prf {
    //! A restraint energy term based on dihedral angles
    /**
      This term provides a way insert an arbitrary energy term that is purely
      a function of dihedral angles. The term can be inserted using the
      settings file or the command line. What angles and what functional
      form to use can be configured with an XML file read in at run time.
      Efficient delta calculations for this energy term are trivial. How to
      use this class to restrain simulations is described in \ref dh_restraints .
      \sa \ref dh_restraints
      */
    class DihedralRestraints : public Energy
    {
    public:
        //! Default constructor
        DihedralRestraints();
        ~DihedralRestraints();
        //! Set parameters from an XML file
        void set_pars(std::string xmlfilename);
        void init();
        double evaluate();
        double deltaE(Update *up);
        void rangeEstimate(double &x1, double &x2);
    private:
        std::multimap<int, RestraintFunction *> c;
        std::string filename;
    };
}
/**
  \page dh_restraints Restraining simulations with a torsional angle potential
  Sometimes we are interested in keeping one of the proteins relatively fixed
  while allowing the other chains to move around freely. Sometimes we have a
  good idea about the secondary structure of a certain region of the chain, and
  want the simulations to not drift away from it too easily, while also not
  making that region entirely immobile. One way to guide the simulations in
  such cases is to use an energy penalty for backbone angles which differ from
  their preferred values. To use such a term in a simulation, you can write
  the following line in the settings file:
  \verbatim
  force_field FF08DR=FF08+Extras:DihedralRestraints(restr_descr.xml)
  \endverbatim
  where \tt restr_descr.xml is an XML file describing the restraints (see below).

  Equivalently, you can specify the use of the restraints in the command line.
  \verbatim
  $ BasicMCRun --add_chain_pdb 1 monster.pdb::A -lcyc -1 -ncyc 100000 \
  --force_field "FF08DR=FF08+Extras:DihedralRestraints(restr_descr.xml)"
  \endverbatim
  Notice that quotation marks are needed around the force_field option in the
  command line.

  The above syntax can be used for any ProFASi simulation program.

  \section specs Configuring the dihedral restraints with an XML file
  A dihedral restraints energy term must be configured so that it knows which
  angles to apply restraints on, and what the restraints look like. This is
  done using an XML configuration file. This is the format of the XML file:

  \li The top level node for this XML configuration file is called
      \tt dihedral_restraints .
  \li Within the scope of the top level node, there should be a series of nodes
  of name \tt restraint .
  \li Each \tt restraint node must have a field called \tt dof_id and a field
  called \tt parameters .
  \li The \tt restraint node can optionally have an attribute \tt type,
  specifying the kind of function to be used for the angle. The default
  function for dihedral restraints is a circular normal form.
  \li The contents of the \tt parameters field depend on the type of function
  chosen with the type field. The parameter names for the default function
  are \tt mean and \tt kappa and weight (N).

  In practice, it would be cumbersome to generate large XML files with such
  a nested structure. But ProFASi's XML module has a template parsing facility
  which can be used to do this much more elegantly. Let's say we have
  generated a long array of "desired" values for the dihedral angles with
  some notion about the uncertainity (width) associated with each angle. We
  have a file like this (let's say the file is called \tt dh_con.table ):
\verbatim
::b:5 2.74191 0.5 10 1 1 2.74191 20
::b:6 -1.81689 0.5 10 1 1 -1.81689 20
::b:7 2.17992 0.5 10 1 1 2.17992 20
::b:8 -1.56731 0.5 10 1 1 -1.56731 20
::b:9 2.62498 0.5 10 1 1 2.62498 20
::b:10 -2.22879 0.5 10 1 1 -2.22879 20
::b:11 1.77849 0.5 10 1 1 1.77849 20
::b:12 -1.39277 0.5 10 1 1 -1.39277 20
::b:13 1.99491 0.5 10 1 1 1.99491 20
::b:14 -1.82772 0.5 10 1 1 -1.82772 20
::b:15 0.678933 0.5 10 1 1 0.678933 20
::b:16 -1.01753 0.5 10 1 1 -1.01753 20
::b:17 2.45219 0.5 10 1 1 2.45219 20
::b:18 -1.22173 0.5 10 1 1 -1.22173 20
::b:19 -0.486947 0.5 10 1 1 -0.486947 20
...
\endverbatim
The format of the above file is entirely up to you. It does not matter which
column you write the DOF identifier string, which column you put the mean
etc. But for the discussion, let's say you want column 7 (counting from 1, like
in awk) to be the mean, column 8 to be the kappa and column 4 to be the weight,
and column 1 to be the DOF identifier string. You then make an XML file with
the following contents (restr_descr.xml):

\verbatim
<dihedral_restraints>
  <formatted_data>
    <format name="restraint">
      <dof_id>$1</dof_id>
      <parameters>
        <mean>$7</mean>
        <kappa>$8</kappa>
        <weight>$4</weight>
      </parameters>
    </format>
    <import_data>dh_con.table</import_data>
  </formatted_data>
</dihedral_restraints>
\endverbatim

ProFASi will see this XML file as follows:
\verbatim
<dihedral_restraints>
  <restraint>
    <dof_id>::b:5</dof_id>
    <parameters>
      <mean>2.74191</mean>
      <kappa>20</kappa>
      <weight>10</weight>
    </parameters>
  </restraint>
  <restraint>
    <dof_id>::b:6</dof_id>
    <parameters>
      <mean>-1.81689</mean>
      <kappa>20</kappa>
      <weight>10</weight>
    </parameters>
  </restraint>
  <restraint>
    <dof_id>::b:7</dof_id>
    <parameters>
      <mean>2.17992</mean>
      <kappa>20</kappa>
      <weight>10</weight>
    </parameters>
  </restraint>
  ...
</dihedral_restraints>
\endverbatim
We did not specify a distribution type, so the von Mises or circular normal
form will be assumed for every restraint. The parameters are as found in the
\tt parameters node.
\sa \ref ds_restraints
  */
#endif // DIHEDRALRESTRAINTS_HH
