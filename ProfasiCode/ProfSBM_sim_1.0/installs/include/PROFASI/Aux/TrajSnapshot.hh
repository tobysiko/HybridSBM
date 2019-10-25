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

#ifndef TRAJSNAPSHOT_HH
#define TRAJSNAPSHOT_HH
#include <vector>
#include <cstddef>

class BinStream;

namespace prf_traj
{
    //! Snapshot of program state
    /**
      This class represents a minimal block of information representing
      the state of an MC run with PROFASI at a certain MC sweep. The snapshot
      contains an identifying unsigned long integer representing the MC sweep,
      an integer index for the temperature, a double precision energy value
      for cross-checks, a couple of numbers representing the state of the
      random number generator and values of all degrees of freedom in the
      system.

      The integer temperature index is unnecessary for a canonnical MC run
      with BasicMCRun. But it is written there anyway to keep a fixed format.
      The index refers to the temperature values in an array stored elsewhere.

      The energy value can be calculated from the coordinates, but storing
      at least the sum in the snapshot gives an easy check that the coordinate
      values are interpreted in the same way when they are read as when they
      are written.

      Two numbers are used to store the state of the random number generators.
      A long integer seed, and an unsigned long integer representing how many
      random number calls had been made until the snapshot point. Using these
      two the random number generator can be restored to the exact same point
      in its history.

      The degrees of freedom values are stored as one single array of double
      precision values. Information about chains in population and the box
      length in the simulation have to be provided separately for the DOF
      info to make sense. The reason why sequence information is not appended
      to every snapshot (to make it self-sufficient) is that it does not change
      during the simulation. A snapshot is a record of the changeables in a
      simulation, which must be augmented with the constants.

      \ingroup prf_trajectory
      \sa prf_traj::Trajectory, prf_traj::TrajSeg

      */
    class TrajSnapshot
    {
    public:
        //! Default constructor
        TrajSnapshot();
        //! Copy constructor
        TrajSnapshot(const TrajSnapshot &t);
        ~TrajSnapshot();
        //! Snapshots can be assigned
        TrajSnapshot &operator=(const TrajSnapshot &t);
        //! Read in properties from a BinStream object
        /**
          The BinStream is assumed to have already read the bytes corresponding
          to the snapshot. In this function, the bytes are converted into
          integers, doubles etc and assigned to the respective fields.
          */
        void fill(BinStream &bst);
        //@{
        /**
          @name Snapshot property access
          */
        inline unsigned long MC_time() const {return icyc;}
        inline double energy() const { return ene;}
        inline int T_index() const {return tindex;}
        inline std::vector<double> & coordinates() { return popconf; }
        inline long random_number_seed() const { return seed; }
        inline unsigned long random_number_ncalls() const { return ncalls; }
        //@}
        //@{
        /**
          @name Setting individual fields
          */
        inline void set_MC_time(unsigned long t) {icyc=t;}
        inline void set_energy(double e) {ene=e;}
        inline void set_T_index(int i) {tindex=i;}
        inline void set_coordinates(const std::vector<double> &crds) {
            popconf=crds;
        }
        inline void set_n_coordinates(size_t i) {popconf.resize(i);}
        inline void set_random_number_seed(long c) {seed=c;}
        inline void set_random_number_ncalls(unsigned long ncl) {ncalls=ncl;}
        //@}
    private:
        int tindex;
        double ene;
        unsigned long icyc,ncalls;
        long seed;
        std::vector<double> popconf;
    };
}
#endif // TRAJSNAPSHOT_HH
