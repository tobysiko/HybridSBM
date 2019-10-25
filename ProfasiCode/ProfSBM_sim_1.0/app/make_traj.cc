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

#include <Aux/Trajectory.hh>
#include <Aux/ProgUtils.hh>

int main(int argc, char *argv[])
{
    prf_utils::ProgArgs optn;
    optn.settings_use_type(0);
    optn.option("output_file","o",1);
    prf::Logger::verbosity=2;
    optn.analyze(argc,argv);
    if (optn.option_given("h") or argc==1) {
        prf::cout<<"Usage:\n\n"
                <<argv[0]<<" -o output_file conf.data_XYZ0 conf.data_XYZ1 "
                <<"conf.data_XYZ2 ...\n\n"
                <<"The program takes a set of binary configuration data files "
                <<"generated in a series of restarts of the same run, and "
                <<"creates a \"traj\" file for the run.\n\n "
                <<"The input files need not be ordered in any particular way. "
                <<"But the paths of the input and output files must be "
                <<"specified relative to the current working directory.\n";
        return 0;
    }

    std::string outfile="output_traj",outdir="";
    if (optn.option_given("o")) outfile=optn.option("o");

    prf_traj::Trajectory traj;
    std::deque<std::string> inputconf;
    for (int i=0;i<optn.n_spare_args();++i) {
        inputconf.push_back(optn.spare_args(i));
    }
    traj.append_list(".",inputconf);
    if (traj.init()) traj.save(outfile);
    else prf::cerr<<"Failed to initialize trajectory with given inputs.\n";
    return 0;
}

/**
  \page traj_gen Generating trajectory files
  \section traj_gen_1 What are trajectory files ?
  PROFASI multi-segment trajectory files are described in more detail in
  \ref traj_files . In short, they are files containing information about
  certain other files. PROFASI stores snapshots of the current state of the
  protein system and random number generator in a binary format during the
  run. The binary data is stored in files called "conf.something". When a
  run is restarted many times, each restart results in a new "conf.specialtag"
  file. The "traj" file is a file which keeps track of these multiple
  "conf.specialtag" files so that information can be querried from the entire
  history of the run together, without actually merging the binary files of
  the segments of the run.

  \section traj_gen_2 What programs generate them ?
  BasicMCRun, SimTempRun, SimAnnealRun, ParTempRun and WLRun generate them
  when they are used for a simulation. The trajectories files are called "traj"
  and are located inside the output directories of the simulation. The program
  make_traj can parse a series of "conf.specialtag" files and generate a "traj"
  file describing them, when the binary files come from a single multi-segment
  run.

  \section traj_gen_2b Analysis of data generated by older PROFASI versions
  If you had a version of PROFASI from 1.1 to 1.4, and have made many
  simulations, the output of the simulations does not have "traj" files. If
  there are "conf.info" files containing layout information about the binary
  conf files, you are on safe grounds. You can run "make_traj" as explained
  below in \ref traj_gen_4 . The program will detect and use the "conf.info"
  files.

  If you have generated data from a PROFASI version so old that it
  did not even generate "conf.info", then things are somewhat more difficult.
  One way is to generate the  missing info files by starting a simulation for
  the same system on the same computer where the original run was made, and
  killing it once it has created the traj file. Along with the traj file, the
  first binary data file will be written with an inline header. It's a binary
  file, but you can open it with "less" and read the header section which is
  human readable text. Copy the text between tags "PROFASI_CONF_HEADER" and
  "END_PROFASI_CONF_HEADER", and paste into a new file called "conf.info".
  That's the file you need to put into the same directory as the binary files
  which you want to analyze.

  \section traj_gen_3 Why does one need traj files ?
  Trajectory files are necessary to extract snapshots from older runs with
  extract_snapshot, reconstruct run-time history with additional measurements
  using extract_props, or even to continue an older run.
  \subsection traj_gen_3a Aren't such things handled by the "conf" files ?
  It was, in PROFASI versions prior to 1.4.8. The meta-data traj file provides
  a cleaner way to handle multi-segment runs. It would be possible to write a
  little code to make sure that extract_snapshot etc continue to work with
  single binary "conf.xyz" files as before, while handling multi-segment runs
  with the traj file as an additional feature. But disabling direct handling
  of the binary files is a more drastic way to draw attention to the preferred
  new way of doing things.

  \section traj_gen_4 Generating a trajectory file from binary conf files
  Do this only if you don't have the "traj" files, i.e., your runs were made
  with a version of PROFASI which did not write traj files during the
  simulations.
  \verbatim
  make_traj -o n0/traj n0/conf.data*
  \endverbatim
  The above takes all files matching n0/conf.data*, parses them to find the
  MC cycle limits, sorts them in the order of the \e starting MC cycles,
  adjust cycle ranges to avoid overlaps, and creates a trajectory file called
  n0/traj. The inputs should be PROFASI binary conf files.

  The binary files should preferably each have their own header sections. If
  they don't, the information of the header section can be provided in a single
  "conf.info" file in the same directory. One "conf.info" file must describe
  the binary files without headers in the same directory.

  The names of the input files (conf...) as well as the output trajectory file
  (n0/traj) should be specified relative to the current working directory.
  make_traj can not process absolute path names.

  \section traj_gen_5 Our reasons for introducing traj files
  In the oldest versions of PROFASI, there was only one file to store system
  snapshots along the run, called "conf" inside the output directories. There
  was no information on the layout of the bytes in the files anywhere. Those
  files can essentially only be interpreted with the same version of PROFASI,
  or if they are re-written with header information.

  PROFASI 1.1 introduced the additional "conf.info" files which contained
  detailed byte-by-byte description of the contents of the binary files. This
  information is now written at the start of the binary "conf" files, in a
  "header" section.

  Another important change was motivated by our experience with running long
  simulations on clusters or supercomputers. When data is transferred
  from a compute cluster to the user's long term storage, the constantly growing
  size of the "conf" files lead to larger and larger quantities of data being
  transferred for later stages of the run. Therefore we decided to start new
  binary snapshot files for every restart of a run, rather than appending to
  the pre-existing "conf" files. With the new split files method, the data to be
  transferred remained roughly constant for all stages of long runs. There were
  now many files now, containing different parts of the information that was
  written to the "conf" files of older PROFASI versions. These files were
  called variously in different development versions of PROFASI, such as
  "conf.bkp_TIMESTAMP", "conf.stage_num" etc.

  The segment-wise storage of configuration snapshots afforded some extra
  flexibility to the user, such as running different parts of the simulation
  on different computers without worrying about how they store the binary
  data. But analysing the data became more cumbersome as scripts now had to
  work through all the snapshot files.  The "traj" file provides a nice way to
  organize all the information in these segment files.

  Another important reason to have the traj file instead of initializing with
  the segment files every time is that it reduces the number of files which
  must be opened and closed during initialisation. The traj file contains
  information on MC cycle limits associated with different segments. So, only
  those segment files are opened which are required.

  \sa \ref traj_files
  */