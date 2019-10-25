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

#ifndef MParTemp_HH
#define MParTemp_HH
#include <mpi.h>
#include "GMC.hh"
#include "../Aux/Permutation.hh"

namespace prf
{
    //! Parallel Tempering
    /**
    * This is the implementation of parallel tempering for ProFASi usingMPI.
    *
    * In this implementation, parallel tempering contains the concept of
    * multiplexing. Given "n" nodes working with "m" temperatures, with
    * n=p*m, with integer p, the system contains m simplexes with p elements
    * in each. At each temperature exchange, a node in a certain simplex
    * attempts to exchange with a node in a neighbouring simplex. But the
    * serial number of the node inside the other simplex is varied randomly at
    * every cycle, so that the system does not segregate into p independent
    * parallel tempering runs with m temperatures each, but evolves together.
    * It is hoped that this leads to better/faster equillibration. The exact
    * extent of its effectiveness is unclear. Of course, if n is set equal
    * to m, it reduces to a simple parallel tempering as before.

    * This class assumes that there is some main program somewhere which has
    * been started with MPI. So, no MPI specific initialisation is done here.
    * But all communication between nodes happens using MPI.
    */

    class ParTemp : public GMC
    {
    public:
        ParTemp(); //!< default constructor
        ~ParTemp();

        //! Enable multiplexing
        inline void enable_multiplexing() {multiplexing=true;}
        //! Disable multiplexing
        inline void disable_multiplexing() {multiplexing=false;}

        //! Attempt a configuration switch, returns 0 if fails.
        /**
         * This attempts an exchange of the temperature between two
         * nodes of neighbouring temperature indices.
         */
        unsigned SwitchTemp();

        int Synchronize(int iatmpts);

        //! Basic setup work
        int Setup();
        int parseCommand(InstructionString s);
        void print_setup();
    private:
        int assign_T_esort();
        int init_T_multiplex_adjacent();
        int init_T_multiplex_staggered();
        size_t nswatmpt,my_rank;
        size_t n_procs,nnodesT;
        std::vector<double> rans;
        bool multiplexing;
        std::vector<size_t> perm;
        unsigned *Tarray;
        int *ndarray;
        double *Erarray;
        
        enum { ADJACENT,STAGGERED,SORTED } initTmethod;
    };
}

#endif

