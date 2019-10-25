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

#ifndef DOF_INFO_HH
#define DOF_INFO_HH
#include <string>

namespace prf {
    typedef enum {
        rigid_body_xyz,
        backbone_torsion_angle,
        sidechain_torsion_angle
    } DOFCategory;
    //! A structure with 8 integers to characterize a degree of freedom
    /**
      A degree of freedom in ProFASi is either a torsional angle, or a
      rigid body coordinate. By convention, information about the overall
      location and orientation of a protein is specified through the
      Cartessian coordinates of the first three atoms on the backbone
      seen from the N-terminus. ProFASi 1.5 provides a unified interface
      to all degrees of freedom in the system. This small class stores
      useful information about a particular DOF in one place, like which
      chain or residue it belongs to. All members are public. It is
      important that these are 8 integers: we intend to build a heavily
      accessed array of a large number of these objects.
      */
    class DOF_Info
    {
    public:
        //! Default constructor
        DOF_Info();
        ~DOF_Info();
        //! Copy constructor
        DOF_Info(const DOF_Info &);
        //! Assignment operator
        DOF_Info & operator=(const DOF_Info &d);

        //! Assign to all fields
        void assign(int kind, int ch_no, int grp_no, int glob_idx,
                    int idx_in_chn, int sp_idx_in_grp, int sp_glob_idx,
                    int sp_idx_in_chn);
        //! Make a string representation
        std::string str() const;
    public:
        //! An enumeration of different kinds of coordinates
        int dof_kind;

        //! Which chain the DOF belongs to
        int chain;
        //! Which group the DOF belongs to (when it makes sense)
        int group;
        //! Global index among all DOFs in the system
        int global_index;
        //! Index among all DOFs in the protein chain
        /**
          This is the answer to the question: what is the index of this
          DOF in this chain ? It's relation to the group_index member below
          is complicated. The protein class may rearrange DOFs from all
          constituent amino acids for convenience.
          */
        int index_in_chain;
        //! Index among DOFs of the same kind in the particular residue
        /**
          This is the index of a particular DOF of a certain kind in a residue.
          There may be rigid body DOFs assigned to a group with group indexes
          from 0 to 8. There could be one or two backbone DOF with group indexes
          starting from 0. This goes for side chain DOF as well. The variable
          specific_index_in_group enumerates DOFs of a certain kind in the
          residue.
          */
        int specific_index_in_group;
        //! Global index among DOFs of the same kind in the population
        int specific_global_index;
        //! Index among DOFs of the same kind in the chain
        int specific_index_in_chain;
    };
}

#endif // DOF_INFO_HH
