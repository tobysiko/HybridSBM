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

#ifndef RandomNumberBase_HH
#define RandomNumberBase_HH
#include <string>

namespace prf
{
    //! Base class for Random number generators for use in PROFASI
    /**
     * It is easy to change the random number generator in PROFASI. Any new
     * random number algorithm can be implemented as a class, and for use with
     * PROFASI, that class has to inherit from RandomNumberBase.
     * \ingroup random_number
     */

    class RandomNumberBase
    {
    public:
        RandomNumberBase();
        virtual ~RandomNumberBase();
        inline std::string Name() const {return thename;}
        inline void Name(std::string stn) {thename=stn;}
        //! Generate a seed from process id, time ...
        int createSeed(int inpseed);
        //! Generates a random number distributed uniformly between 0 and 1
        virtual double shoot();
        //! As above, but between a and b
        inline double shoot(double a, double b) {return a+(b-a)*shoot();}

        //! Between 0 and width
        inline double shoot(double width) {return width*shoot();}

        //! Create a random bit, 0 or 1 with equal probabilities
        inline int shootBit() {return ((shoot()<0.5)?0:1);}

        virtual int ConfSize();
        virtual std::string ConfSignature();
        //! Write to a file
        virtual void WriteTo(FILE *cnfil);
        //! Read in the state of the generator from a file
        virtual void ReadFrom(FILE *cnfil);
        //! Write to a file in text form
        virtual void WriteTo_text(FILE *cnfil);
        //! Read in the state of the generator from a file in text form
        virtual void ReadFrom_text(FILE *cnfil);
        //! Randomize State based on a given "rank".
        /**
         * This could be used to make sure that different ranks in a parallel run
         * start with different random number states.
         */
        virtual void RandomizeState(int irnk=0);
        //! Initialize with a given seed
        virtual void ResetDefaultState(long seedd);
        //! Retrace history
        /**
        * Each random number object has a "seed" from which it starts
        * calculating. Every time a random number is generated, a counter is
        * incremented, so that at all times, the number of calls to the random
        * number generator since the beginning is known. These two numbers can be
        * used to restore the random number generator to a particular state. This
        * function causes the random number generator to first set its seed to
        * what it knew for a seed, and then retrace all its history by
        * regenerating an equal number of random numbers. With this function, one
        * can save only two numbers in the configuration files for a random number
        * generator, and restore the state as required. This restoration procedure
        * is somewhat slower, but it is much better from the point of view of
        * disk-space in configuration files. So, from PROFASI version 1.1, this
        * will be the default way of saving and restoring states of a random
        * number generator.
        */
        virtual void retrace(std::string sttfile="");
        //! Set seed and num_calls (to later call retrace)
        inline void set_seed_and_ncalls(long sd, unsigned long ncl) {
            myseed=sd;
            ncalls=ncl;
        }

        //! Retrieve number of calls
        inline unsigned long num_calls() {return ncalls;}

        //! Sets seed and Initializes with it
        virtual void setSeed(long i);
        //! Returns the seed
        /**
         * Note that this is not necessarily the "current seed". For some
         * random number generators, it may be possible to capture the entire
         * state conveniently with a seed, but it is not assumed here for all
         * random number generators. The seed here only makes sure that
         * different seeds mean different sets of generated values.
         */
        inline long Seed() const {return myseed;}

        //! Save complete state for faster recovery
        /**
          * As runs get longer and longer, one can not afford to lose the
          * first 15 minutes to restoring the random number generators to
          * the state at the end of the previous run. Especially if the
          * runs are limited to 6 or 12 hours each time, on some supercomputer.
          *
          * Therefore, as a compromise between disk usage and restart speed,
          * a random number generator class in PROFASI will, from version 1.4,
          * implement a saveState() function. This function will store in a
          * text file all internal state variables of the generator. This
          * function is meant to be called at the end of each small run, so that
          * the disk usage is not a concern. But this will ensure instant
          * restart, if the run is to resume from the end of a previous run.
          *
          * The format of data in the state files must contain the name of the
          * random number generator, the regular information like the initial
          * seed and the number of calls. The content of the files after this
          * point could be specific to the random number generator class.
          */
        virtual void saveState(std::string sttfile);
        //! Recover complete state from file
        /**
          * Check if the given file has a state stored by a random number generator
          * object with the same name. Check if the original seed written in the
          * file matches with what is already stored. Check if the saved state in the
          * file is in the past relative to the number of calls. If all checks pass,
          * read in the state from the file, and if necessary, retrace the remaining
          * number of random number generations to fast forward the ncalls to the
          * required number.
          */
        virtual int recoverState(std::string sttfile);

    protected:
        long myseed;
        unsigned long ncalls;
        std::string thename;
    };
}

#endif
