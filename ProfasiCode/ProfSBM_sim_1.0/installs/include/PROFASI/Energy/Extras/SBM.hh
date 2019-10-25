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

#ifndef SBM_HH
#define SBM_HH
#include "../Energy.hh"
#include "RestraintFunction.hh"

namespace prf {
    //! A single distance restraint
    /**
      A distance restraint consists of a pair of atoms and a function
      of the distance between them.
      */
    class SBMrestraint {
    public:
    	SBMrestraint();
        ~SBMrestraint();
        void delete_restraint_function();
        SBMrestraint(const SBMrestraint &dr);
        SBMrestraint &operator=(const SBMrestraint &dr);
        int set_pars(prf_xml::XML_Node *pars, Population *p);
        double evaluate();
        double estimate_upper_bound();
        double estimate_lower_bound();
        inline void setPerRestraintEnergy(double e){perRestraintEnergy=e;f->setPerRestraintEnergy(e);}
        inline double getPerRestraintEnergy(){return f->getPerRestraintEnergy();}
        inline bool is_in_range(int first, int last){ 
        	if ( res1>=first && res1<=last && res2>=first && res2<=last ) 
        		return true; 
        	else 
        		return false;
        }
        inline int getRes1() {return res1;}
        inline int getRes2() {return res2;}
        
    private:
        int get_aid(std::string atm, Population *p, int &r);
        int atom1,atom2,res1,res2;
        RestraintFunction *f;
        double perRestraintEnergy;
       
        
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
    class SBM : public Energy
    {
    public:
        //! Default constructor
    	SBM();
    	SBM(std::string m);
    	SBM(double ne);
        ~SBM();
        //! Set parameters from an XML file
        void set_pars(std::string xmlfilename);
        void init();
        double evaluate();
        double deltaE(Update *updt);
        void rangeEstimate(double &x1,double &x2);
        inline void setNormEnergy(double e){normEnergy=e;}
        inline double getPerRestraintEnergy(){return perRestraintEnergy;}
        inline int getNrestraints(){return nRestraints;};
        inline void set_coop(double val){coop = val;};
        inline double get_coop(){return coop;};
        /*bool isContactInRigidRange( int cind, int s, std::vector<std::pair<int,int> > *rng, Update *updt);*/
        inline void setDebug(bool b){debug = b;}
        inline void setScaleSBM(double s){scale_SBM = s;}
        inline double getScaleSBM(){return scale_SBM;}
    private:
        std::vector<SBMrestraint> c;
        std::string filename, posmapfilename;
        std::string restraintID;
        int nRestraints;
        double normEnergy;
        double perRestraintEnergy;
        double coop;
        int nchanges;
        bool debug;
        double scale_SBM;
        
    };
    
    
    
}
/**
  \page sbm_restraints Simulations with SBM (Structure Based Model) native contact distances
Using SBM restraints in a simulation follows exactly the same syntax
as does the use of distance restraints (See \ref ds_restraints ). The following
one line in the settings file will set up a new energy term according to the
restraints specified in the XML configuration file \tt dist_restr.xml .
\verbatim
force_field FF08DR=FF08+Extras:SBM(sbm_restr.xml)
\endverbatim
The procedure for specifying the restraints in the command line is the same as
for the distance restraints.

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
function for SBM is FMULTIGAUSS, i.e. a fixed-depth multi-well Gaussian.
\li The contents of the \tt parameters field depend on the type of function
chosen with the type attribute. The parameter names for the default function
are \tt mean and weight.

\section example Example
The following XML file configures 5 distance restraints:
\verbatim
<sbm rid="_superGA_prf_min_rmsd">
<formatted_data>
  <format name="restraint" type="$3">
  <atom1>$1</atom1>
  <atom2>$2</atom2>
  <parameters>
    <minima>$4</minima>
    <radius>$5</radius>
    <steepness>$6</steepness>
    <width>$7</width>
    <depth>$8</depth>
  </parameters>
  </format>
  <data>
    0/0/MET/_CA_  0/4/ASP/_CA_  FMULTIGAUSS  7.839240652,7.73980858936,7.49781708232,6.41458852928  0.0  1.0 0.5 1.0
    0/7/SER/_CA_  0/11/ALA/_CA_  FMULTIGAUSS  6.63486103848,5.98854523236,7.08288874683,8.16002518378  0.0  1.0 0.5 1.0
    0/7/SER/_CA_  0/12/LYS/_CA_  FMULTIGAUSS  9.15407920001,8.58863702807,9.19692932451,9.36218852619  0.0  1.0 0.5 1.0
    0/7/SER/_CA_  0/38/VAL/_CA_  FMULTIGAUSS  7.99676834728,6.53995856256,7.0349295661,7.89524705123  0.0  1.0 0.5 1.0
    0/8/LEU/_CA_  0/12/LYS/_CA_  FMULTIGAUSS  6.85261278638,5.75008773846,6.0492036666,6.31429093406  0.0  1.0 0.5 1.0
  </data>
</formatted_data>
</sbm>
\endverbatim

In connection with distance restraints, we introduced the template interpretation
capacity of the ProFASi XML module. The above is another example. Here the
tabular data to be re-cast as a series of XML nodes appears inline inside a
special tag \tt data. The \tt format tag is used to interpret the contents
of the data tag. Here we also introduce how to specify the functional form
by giving a type attribute to the auto-generated restraint nodes.

The atoms are specified through their chain, residue number, residue type, and
atom label. The labels are 4 characters wide. Spaces in the label should be
replaced by underscores. It is important to have the labels correct. It is
" CA " (and therefore "_CA_" above), and not "CA__" for C-alpha.

Several alternative forms are available for the functional form. [TODO]

\sa \ref ds_restraints
  */
#endif // SBM_HH
