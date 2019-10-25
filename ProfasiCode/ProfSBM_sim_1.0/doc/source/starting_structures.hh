#ifndef STARTING_STRUCTURES_HH
#define STARTING_STRUCTURES_HH

/**
  \page initstruct Initial structures in PROFASI simulations
  Almost all our simulations start with random initial values assigned
  to all degrees of freedom. In a few instances we have used a "stretched"
  chain as the starting point. It is however a frequent request that it should
  be possible to start simulations from a given initial structure. Here we
  describe how that can be done.

  \section initstruct_ranstr Random and stretched initial conformations
  To start simulations from random configuration of the protein chains, you
  don't have to do anything about the start configuration. That is the default
  behaviour.

  To start from stretched chains, use the command \e init_config in the
  settings file.
  \verbatim
  add_chain 1 <*GEWTY DDATKT FTVTE*>
  init_config stretched
  \endverbatim

  Simulation programs based on BasicMCRun, such as SimTempRun, SimAnnealRun,
  ParTempRun and WLRun store the starting structure for the simulations in a
  file called \e start.xml. Each rank in a parallel run stores its own start
  conformation. Normally, the start conformations are stored after the so
  called "preliminary relaxation" cycles. These are normal MC cycles at a
  high temperature to remove atom clashes, so that in the main loop of the
  simulation, the energies start at reasonable values. The number of relaxation
  cycles is 10 times the number of chains in the system. This is the reason
  why even when you specify "stretched" as init_config, you don't see a flat
  elongated chain as the start structure.

  Remember that "stretched" does not really mean using a physical stretching
  force and determining the stretched protein conformation. It simply assigns
  values to the degrees of freedom such that the chain appears extended. It
  is not guaranteed that a stretched conformation in this sense wont lead to
  clashes between atoms. Therefore the relaxation cycles are run even for
  stretched starts. It is possible to use a settings file command to skip the
  relaxation cycles, so that init_config stretched will produce a stretched
  start.xml file.

  \verbatim
  add_chain 1 <*GEWTY DDATKT FTVTE*>
  init_config stretched
  preliminary_relaxation off
  \endverbatim

  \section initstruct_given Starting from a given structure
  Here we are talking about starting the simulations from, for instance, a PDB
  file. Suppose you have a PDB file \e abc.pdb and you want a simulation to
  start from that structure. Follow these steps:

  \li Generate a regularized approximation of the structure with PROFASI's
  constraints, like this:
  \verbatim
  regularize abc.pdb
  \endverbatim
  This will produce two files "min_etot.xml" and "min_rmsd.xml". Check how good
  the approximation is by converting them to PDB files :
  \verbatim
  prf_convert min_etot.xml min_etot.pdb
  \endverbatim
  Use your favourite molecular visualisation program to superimpose min_etot.pdb
  with abc.pdb and convince yourself that a good approximation has been found.
  \li Rename either the min_etot.xml or min_rmsd.xml to something convenient
  like "abc_init.xml"
  \li Scrap the add_chain or add_chain_pdb command in the settings file and
  do this in the settings file instead:
  \verbatim
  set_population abc_init.xml
  \endverbatim
  The "preliminary relaxation" cycles are not run if there is an explicitly
  mentioned initial structure.

  The above illustrates the principle. The recommended way to provide a starting
  structure is to use a structure in \ref xmlstruc. The programs
  \e regularize and \e prf_convert can be used to generate such a structure from
  a given PDB file. Regularize does more than just the format conversion, and is
  recommended unless the PDB file was generated in a PROFASI simulation. The
  XML structure file can also be extracted from a previous simulation with
  PROFASI using the program \e extract_snapshot .

  \sa \ref regul , extract_snapshot , \ref prf_convert, \ref xmlstruc
  */
#endif // STARTING_STRUCTURES_HH
