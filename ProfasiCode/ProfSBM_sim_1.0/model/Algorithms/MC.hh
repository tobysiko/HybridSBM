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

#ifndef MC_HH
#define MC_HH
#include <sstream>
#include "../Aux/fileutils.hh"
#include "../Updates/UpdatesHandler.hh"
#include "../Energy/FFHandler.hh"
#include "../Elements/Population.hh"
#include "../Aux/Constants.hh"
#include "../Aux/HandlerBase.hh"

/**
  * \defgroup MCalgorithms Monte Carlo algorithms
  * Markov Chain Monte Carlo algorithms
  *
  * ProFASi provides several Monte Carlo algorithms such a simulated tempering,
  * parallel tempering, simulated annealing. They are based on a class called
  * MC, for Monte Carlo, which implements elementary Metropolis-Hastings
  * Markov Chain generation. But MC can also be sub-classed to implement
  * many algorithms such as the Multicanonical sampling or the Wang-Landau
  * algorithm which have a different accept/reject condition. The classes
  * in this module only provide the background. To see how to do this kind
  * of simulations in ProFASi, see the application programs BasicMCRun,
  * SimTempRun, SimAnnealRun or ParTempRun. The tutorials also contain a
  * basic description of how to perform such simulations.
  *
  */
namespace prf
{
    //! MC class handles the Markov Chain Monte Carlo evolution
    /**
     * The MC class accumulates information about the force field and updates
     * and the population of proteins. It then performs the Markov Chain
     * evolution of that system, using predefined probabilities for the
     * updates. MC as such performs only fixed-temperature simulations. For
     * most practical purposes this would be inadequate. But it works as a good
     * base class for other methods such as the Simulated Tempering (SimTemp
     * class) or Simulated Annealing (SimAnneal class).
     *
     * The acceptance probability for a proposed update is very much like
     * the Metropolis probability, but with a prefactor. For certain updates,
     * like BGS, it is necessary to use an additional factor in the probability
     * to ensure detail balance.
     *
     * In version 1.5 handling of energy terms is delegated to a new force field
     * handler class, FFHandler. Similarly operations relating to the conformational
     * updates are delegated to the UpdatesHandler. Both these handlers are
     * members of the MC class. Using the handler classes help maintain a modular
     * structure in the code. These measures are also the first step towards
     * making MC algorithms such as parallel tempering, simulated tempering or
     * simulated annealing available as a pure template library, independent of
     * the rest of ProFASi. Such a library could be useful for arbitrary systems
     * not necessarily related to proteins.
     *
     * Another change relative to earlier versions of ProFASI, made to bring the
     * MC class closer to a pure representation of the Markov chain Monte Carlo
     * algorithm was to rid it of the responsibility of keeping track of the
     * number of temperatures and the temperature index for simulated tempering,
     * parallel tempering and simulated annealing. Although these variables are
     * common to these three algorithm, they do not belong to the concept of
     * the basic Monte Carlo algorithm.
     * \ingroup MCalgorithms
     */

    class MC : public HandlerBase
    {
    public:
        MC(); //!< default constructor
        virtual ~MC();

        //! Return pointer to the manager of conformational updates
        inline UpdatesHandler *updates_hander() { return &uph; }

        //! Return pointer to the force field manager
        inline FFHandler *forcefield_handler() { return ffh; }

        //! Introduce the force field manager
        inline void forcefield_handler(FFHandler *fh) { ffh=fh; }

        //! Return a pointer to the random number generator for other use
        inline RandomNumberBase *RandomNumberGenerator() { return ran; }

        //! Force use of a different a random number generator
        inline void RandomNumberGenerator(RandomNumberBase *grn) {ran=grn;}

        //! Get temperature in dimensionless ProFASi units
        inline double temperature() const {return 1.0/bt;}

        //! Set temperature in dimensionless ProFASi units
        inline void temperature(double x) {bt=1.0/x;}

        //! Get MC temperature in Kelvin
        inline double temperature_in_kelvin() const {return UnivConstants::pru_in_kelvin/bt;}

        //! Assign a temperature in kelvin for the simulation
        inline void temperature_in_kelvin(double val) {
            SetBeta(UnivConstants::pru_in_kelvin/val);
        }

        //! Set inverse temperature (model units for temperature)
        void SetBeta(double be);

        //! Tell MC what population of proteins it needs to evolve
        void Connect(Population *p);

        void set_output_prefix(std::string prfx) {output_prefix=prfx;}

        //! Initialize
        /**
          * After this function, the MC object should be ready to do Monte Carlo
          * steps.
          */
        virtual int Setup();

        //! Perform one Monte Carlo step: naive version
        /**
         * Select an update. Find energy change. Ask the Metropolis question.
         * Decide whether or not to accept the update, and do so.This function,
         * is not used in ProFASi, except during debugging.
         */
        virtual int SimpleStep();

        //! One Monte Carlo Step, as used in PROFASI
        /**
         * Similar to the SimpleStep function, but with a speed up trick. The
         * random number to be used for the Metropolis question is completely
         * independent of the energy calculations. In this version of the
         * Monte Carlo step, that random number is generated before the energy
         * calculations (but after the conformational update). Using it, and
         * the temperature, a maximum acceptable energy increase is calculated
         * before doing anything with the energy terms. This value is passed to
         * the force field handler for the delta E calculations. It is
         * understood that if the energy change is greater than this limit,
         * it need not be calculated precisely. The conformational update will
         * be rejected anyway! So, the FFHandler gets a license to return any
         * (too) large value, if it has ascertained that the correct value lies
         * beyond the limit passed to it. See the documentation of the
         * deltaE_with_limit() function in FFHandler for a little more on this.
         */
        virtual int Step();

        //! Perform nstps Monte Carlo steps
        inline void Step(int nstps) {for (int i=0;i<nstps;++i) {Step();if (debug) std::cout<<"\nSTEP "<<i<<"\n";}}

        //! Run a Monte Carlo cycle or sweep
        /**
         * A Monte Carlo sweep or cycle consists of a pre-configured number of
         * Monte Carlo steps. The default number of steps ( ProFASi v 1.5 )
         * is the number of degrees of freedom in the population. This may be
         * changed to any value using the command "steps_per_cycle N". This
         * function is a short-hand to call Step() many times.
         */
        inline void RunCycle() {
        	if (debug) std::cout<<"MC.hh  inline RunCycle(); nSteps="<<cyclgt<<"\n";
            Step(cyclgt);
            popl->EnforceBC();
            popl->Reconstruct();
            AtomCoordinates::update(0,popl->NumberOfAtoms());
        }

        //! Set length of cycle to i
        /**
         * If a positive value is passed, cycle length is set to that value.
         * If a negative number or 0 is passed, cycle length is set equal to
         * the number of degrees of freedom in the system.
         */
        void SetCycleLength(int i);
        inline int CycleLength() const {return cyclgt;}
        void setup_cycles();

        //! Execute an instruction
        /**
          * The MC class handles some instructions of its own, as well as
          * forwarding instructions meant for the UpdatesHandler. Instructions
          * for the MC class are:
          * \li <b> steps_per_cycle : </b> Sets how many Monte Carlo updates
          * constitute a Monte Carlo cycle. Example: steps_per_cylce 100
          * \li <b> random_number_seed : </b> Manually set a random number seed.
          * \li <b> temperature </b> Temperature for the MC simulation. Example:
          * temperature 310 Kelvin, or temperature 0.485. If no unit is given,
          * temperature is read in the dimensionless ProFASi units.
          * \li <b> ntmp </b> Number of temperatures (mostly relevant for
          * methods derived from Metropolis Monte Carlo) Example: ntmp 10
          */
        virtual int parseCommand(InstructionString s);
        virtual void parseCommands(std::list<InstructionString> &cmds, int argc, char *argv[]);
        virtual void show_help();
        virtual void print_setup();
        virtual std::string ConfSignature();
        inline void setDebug(bool b){debug = b;}
    protected:
        virtual Update * perform_update();
        UpdatesHandler uph;
        FFHandler *ffh;
        Population *popl;
        int cyclgt;
        std::string output_prefix,descr;
        RandomNumberBase * ran;
        double bt;
        bool debug;
    };

}

#endif
