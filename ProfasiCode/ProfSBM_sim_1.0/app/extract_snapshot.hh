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

#include <fstream>
#include <sstream>
#include <Elements/PopulationHandler.hh>
#include <Aux/ProgUtils.hh>
#include <Aux/Trajectory.hh>

using namespace prf;
using namespace prf_utils;

//! Extracting a snapshot configuration
/**
 The program extract_snapshot browses through the saved trajectory information
 from a run of a typical profasi simulation program like BasicMCRun,
 SimTempRun or ParTempRun, and extracts the state of the population at a
 requested  point. The output could be a pdb file, an xml PROFASI structure
 description file, a text or binary structure file with only the population
 information.

 Typically, because of job time limits on clusters and supercomputers, a
 simulation requires many restarts. Every stage of such a run opens a new
 binary data file with a name "conf..." inside the run directoy (n0,n1 etc.).
 In addition, with PROFASI version 1.5, there are trajectory metadata files
 called "traj" in the run directories. They store information about ranges,
 write intervals and number of blocks present in different conf files. For
 extraction of structures and properties, it is best to use these "traj" files,
 as input to this program. The traj files contain information about which
 binary file to open to look for the desired structure. But since they don't
 have any data as such, if you move the traj files, you should move all the
 corresponding "conf..." files along with it.


 \note PROFASI's binary configuration files store information
 in an assumed order. So, this program will not work if the configuration
 layout is changed. The original configuration files typically include random
 number generator information etc. in addition to the state of the population.


 If the reason for extracting a structure is to start new runs using it,
 we recommend using the xml structure file for this, as it contains sequence
 information, and works as a self contained human readable record of the
 structure, unlike the textconf or the binary configuration files. Runs
 can be started with the xml configuration by using "set_population somefile.xml"
 command, while with the text configuration files, you would need to use
 the "add_chain ..." commands as well as "init_config" command. (See \ref
 settings_population for more information on these commands.) \n \n

\section syntax Syntax

extract_snapshot [OPTIONS] input_trajectory_file \n \n
Options could be any of ...\n
    \li <b>-o  or --output_file </b> 1 argument, the filename
    \li <b>-f or --output_format </b> Here you choose whether the output is
    written in pdb, xml, text or  binary configuration format.
    Possible values are pdb, xml, binary and textconf.
    \li <b>-c or --cycle_number </b> The Monte Carlo cycle number at which
   extraction is desired.
    \li <b>-i or --stdin </b> For extracting at given list of cycle numbers. In
 this mode, the program waits for the input of a list of cycle numbers from
 standard input, and when it reads EOF, it extracts a list of snapshots from
 each of those points. The cycle numbers are appended to the end of the
 corresponding output files. The intended use is to "pipe" in a set of cycle
 numbers to the extract program.
    \li <b>-r or --raw </b> This option is provided so that it is possible to
    work without the metadata file "traj", directly using the data. But it is
    inefficient and should be avoided.
\section examples Examples

\verbatim
$ extract_snapshot -o min_at_1293999.pdb n19/traj -c 1293999
$ extract_snapshot -o min_at_1293999.xml n19/traj -c 1293999
$ extract_snapshot -o min_at_1293999.xml -c 1293999 --raw n0/conf.bkp0 conf.bkp2 conf.bkp1
\endverbatim

In the above examples, we extract the state of the population at the cycle
1293999 and save the state as a pdb or xml file. The last example shows that
it is possible to use the binary conf files directly with the raw option. It
is possible, but not recommended for efficiency. The format of the data in the
output file can be inferred from the suffix you provide. So, the "-f" option
is unnecessary. In the following case, the "-f" option would be necessary:

\verbatim
$ extract_snapshot -o weird_structure.whatonearth n0/traj -c 87999 -f pdb
\endverbatim

Without the "-f pdb" option above, the file "weird_structure.whatonearth" would
be written in the (default) xml format. When given, the option "-f" overrides
any automatically inferred output file format.

In order to interpret the binary data in the conf files, information about
the layout of data in these files is needed. From ProFASi version 1.5, this
information is embedded in the binary files. In versions 1.1 to 1.5, this
information was in a separate file called "conf.info" in the same directory
as the conf file. For hints on extracting information from old runs, made with
PROFASI versions without embedded layout information or conf.info files, see
the section "Analysis of data generated by older PROFASI versions" in
\ref traj_gen.

The following (BASH) shell script will create files emin.xml corresponding
to the minimum energy configurations from each of the runs n0...nN, and save
them in the respective directories. Those population configurations can then
be used to start a run with different initial configurations for each of the
nodes. \n \n

\verbatim
$ for i in `seq 0 N` ; do
> tx=`grep ENERGY n\$i/minen.pdb | awk '{print \$9;}' ` ;
> extract_snapshot n\$i/traj -o n\$i/emin.xml -c \$tx ;
> done
$
\endverbatim

This works because PROFASI saves energy and Monte Carlo time information as
remarks in its pdb output files when it can. Note that some of the quote marks
above (used around the seq and grep commands) are back quotes. <br><br>

Finally, note that if you see this line in the documentation, your version of
extract_snapshot is capable of reading binary configurations written in other
machines with different binary formats. If you made a run on a cluster of IBM
PPC 6 processors, and want to analyze the data on your own laptop
with an old Intel Pentium M in it, you don't have to convert the binaries.
Use it directly:

\$ extract_snapshot -o min_at_1293999.xml --raw n19/conf -c 1293999 \n \n

The byte order information is part one of the things written in the header part
of the binary files in recent versions and in the conf.info file in older
versions of PROFASI. In the above, using the raw option is not inefficient,
as there is only one conf segment. If the cycle 1293999 is found in it, it
will be recovered.

\sa extract_props
*/

class extract_snapshot
{
public:
    extract_snapshot();
    ~extract_snapshot();
    int my_init(int argc, char *argv[]);
    int execute();
private:
    void show_help();
    int act_on_CL_options();
    int add_segments(std::string seginfo);
    PopulationHandler PH;
    ProgArgs optn;
    std::string outfile, oformat, settingsfile;
    std::vector<std::string> formats;
    unsigned long cycno;
    prf_traj::Trajectory traj;
};
