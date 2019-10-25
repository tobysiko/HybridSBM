/**
* \page tutorial_1 Tutorial 1: A half an hour tour of ProFASi
\section prereq What you need for this tutorial
If you have just downloaded ProFASi, successfully compiled it, and are wondering
what it can do for you, this is the right place to begin. Even if you have had
the program for some time, you can get a tip or two from this page. If you have
not yet been able to compile the program, you should look at the
\ref installation .

You should be reasonably familiar with the UNIX command line to go through the
following exercises. ProFASi has no graphical user interface. To visualize the
structures generated by the program, you will need a molecular visualization
program, such as <a href="http://pymol.sourceforge.net/">PyMOL</a>.

To start with, create a fresh directory for the following tutorial and
change directory into it. In the following, we will refer to this directory
as your "run directory", and the location of your ProFASi applications directory
as \$profasi_dir. If you set up an environment variable profasi_dir to point to
the subdirectory app under your ProFASi directory, you can copy and paste commands
from the tutorials into your shell. For instance, consider a user "crocodile",
who uses the BASH shell, and has downloaded ProFASi to the directory
/home/crocodile/profasi, and successfully compiled the program. For the
tutorials, it would be sufficient to write, \n

\c crocodile\@river:/home/crocodile> export profasi_dir=/home/crocodile/profasi/app

\section basicmc A simple Monte Carlo simulation
Let's start by running a simple Monte Carlo simulation on a 37 residue peptide
to see how things work. Type the following command in your run directory:

\verbatim
prompt> $profasi_dir/bin/BasicMCRun \
-ac 1 "* DTASDAAAAAALTAANAKAAAELTAANAAAAAAATAR * NH2" \
-T "300 Kelvin" -ncyc 10000 -nrt 1000 -iavg 1000
\endverbatim
This should start a canonical Monte Carlo (MC) simulation with one chain of the
peptide with sequence "DTASDAAAAAALTAANAKAAAELTAANAAAAAAATAR-NH2" at temperature
300 Kelvin. It will perform 10,000 MC cycles and write a record of various system
properties into a run-time history file (normally called "rt" by ProFASi simulation
programs) at an interval of 1000 cycles. It will also calculate averages of
various measured quantities every 1000 sweeps. Let it run. It wont take very long!

\subsection sweep_def Monte Carlo steps vs sweeps
At this point, we would like to clearly define what we mean by a sweep or a Monte
Carlo cycle as opposed to a Monte Carlo step. An MC step consists of a
conformational update being proposed and accepted or rejected based on a
Metropolis criterion. A sweep or cycle consists of a certain fixed number of
elementary Monte Carlo updates or steps.

The default size of a sweep in ProFASi 1.5 is equal to the number of degrees of
freedom in the system. But the  option "-lcyc" (length of a cycle) can be used to
set it to any other value.

\subsection cmd_opts Command line options
The option "-T" is a short form for --temperature and "-ac" means "--add_chains".
It is necessary to explicitly write "Kelvin" as above. Just typing -T 300 would
set temperature to 300 units in PROFASI's internal unit for temperature. 666.67
Kelvin corresponds to 1 in this internal temperature unit. So, 300 units would
be a somewhat high temperature to simulate protein folding.

The first argument to "ac" is the number of chains and the second is the amino
acid sequence. You should look in the section \ref seq_input to learn more about
how to write peptide sequences for different ProFASi programs.

Across all ProFASi programs that accept options, the options can be specified in
an arbitrary order, and not all of them must always be given. So, the user does
not need to worry about what order to specify options on the command line. Also,
most often one can obtain a list of available options by typing the name of the
program without any arguments, or running it with only a "--help" option.

\subsection output_files Output generated by the example simulation
So, what happened in the simulation ? You will notice that the program created a
new sub-directory called "n0". All output files produced by the run are in that
directory. We will now describe them.

First the name of the sub-directory itself. It is called n0 because it is a
single process. When you run the ProFASi simulation programs in the parallel
mode, the output files from different computation nodes are saved in
sub-directories n0,n1,n2 ...

\subsection pdbfls Current and minimum energy structures
The file "minen.pdb" (in n0) is the minimum energy conformation found in this
run. You can view it in your favourite PDB viewer. Even if your run lasted much
less than a minute, take a look! The correct native structure for this small protein
is a single alpha-helix. It is very likely that your run yielded a structure with
some helix content. "current.pdb" is the last structure saved in the simulation.
The structure just before the main Monte Carlo evolution starts (but after a few
high temperature relaxation cycles have been run to get rid of steric clashes) is
stored in the file "start.xml" (which can be converted into a pdb file using the
program prf_convert. See \ref prf_convert ).

\subsection optrs Random numbers
Now, remove the directory "n0" and run the same command again, and look at the
minimum energy and last structures. You will see that every time you run the
program, it produces different minimum structures. This is because by default
these programs start with random initial conformations with a random number seed
derived from the UNIX process id for that run and the current time. Practically,
every run is independent of every other. If you want to force the program to
trace the same trajectory, you have to give it the same random number seed. In
every run, the program writes the random number seed to the standard output. Look
for it from the output lines printed during the example run and run the command
again with an additional option "--random_number_seed that_seed". The program will
retrace its path and produce the same minimum energy conformation.

In a short run like this, you will not always find that the minimum energy
structure ("minen.pdb") resembling the native structure of the protein.  But if
you run it many times, you are likely to find helices more often, and in a
long simulation, probably lasting a day or so, you will be able to determine the
statistical weights of different kinds of structures. It is our educated guess
that you did not expect publication worthy statistics from a sub-minute
simulation.

\subsection timelim Time limited runs
Most serious folding simulations are done in some form of high performance
computing facility, where a simulation is not allowed to run indefinitely. There
are time limits for the jobs, and if a job exceeds this limit, it is killed by
the batch system. ProFASi 1.5 simulation programs are aware of the available
time. Clear the output directory n0 and rerun the simulation as follows:

\verbatim
prompt> $profasi_dir/bin/BasicMCRun \
-ac 1 "* DTASDAAAAAALTAANAKAAAELTAANAAAAAAATAR * NH2" \
-T "300 Kelvin" -ncyc 10000 -nrt 1000 -iavg 1000 -time 00:01:00
\endverbatim
The "time" option tells the program that there is a time constraint, and if it is
not through with its task by then, it will suspend cleanly just before that time
limit expires. When it suspends, the program writes a message to the standard
output about how long it thinks it will need to finish the job. Suppose, on
your computer it says that it needs 3 more minutes. How do you continue the run
from where it suspended ? Easy. Just type the above command again. It detects
that there was a previous incomplete simulation, and it continues from the
end of the previous simulation.

\subsection cctn Caution
It is important to remember that a canonical Monte Carlo simulation, like the one
performed by the BasicMCRun program, can easily get stuck in a local minimum of
the energy landscape. The computation time required to properly sample the energy
landscape of even a small peptide like this example is unnecessarily large. That
is why there are more sophisticated Monte Carlo techniques. Some of them are
implemented and provided as application programs in ProFASi. The goal of this
tutorial is not to get a correct folding event, but to get familiar with different
compoents of ProFASi.

\subsection extrct Extracting structures and the program configuration file "conf"
Ok! Perhaps your simulation resulted in a minimum energy structure with a well
formed helix, or perhaps not.  But how did it get there ? ProFASi simulation
programs do not save a PDB snapshot after every sweep. If we did, you would have
ended up with 10000 PDB files in a minute with the above run!  But it is
desirable to be able to generate the PDB snapshot at any chosen point in a
finished simulation during a post run analysis. To help with this, the simulation
programs like BasicMCRun store their state to the disk at regular intervals.
This interval is normally 1000 Monte Carlo <a href="#sweep_def">sweeps</a>, but
can be altered with the option "-iconf conf_write_interval" to BasicMCRun.

In PROFASI 1.5, the data about the program state is written in several related
files. The topic is discussed in detail in \ref traj_files. Briefly, there is an
"envelope" file called \e traj in your output directory (n0 here), which contains
information about other files where trajectory data is written by the same run.
These other files are a series of binary "conf.dataXXXX" files, where the "XXXX"
stands for a timestamp for the first write into that file. If a run has to be
done in 53 time-limited segments of 6 hours each (due to job time limits on the
cluster or supercomputer) each of the 53 segments will start a new "conf.dataXXXX"
file with the appropriate time stamp in its name, and register it with the
envelope traj file. This way, it is possible to read the (plain text) traj file
and tell which conf.dataXXXX file ought to be opened to get information about a
certain MC cycle. The analysis programs provided with PROFASI 1.5 do this for you
so that you don't need to worry about the mechanism. But you should understand
that although it will seem like you use the "traj" files to extract information,
traj file is only like a telephone directory. The real data is in the much larger
binary conf files. So, if you decide to copy the results of your run, you have to
copy all the files in your n0, n1 ... directories, not just the traj files.

There is a helper application to extract the protein states from a run:
extract_snapshot. Type \n
\c prompt\> \$profasi_dir/bin/extract_snapshot \n
in the present directory. You will be presented with some usage information. This
is a very useful application, so, please spend a minute getting familiar with it.
Extract the PDB file corresponding to the different MC-times during the
simulation at which the program state was saved. If iconf is set to 1000 (default),
the configurations will be saved at cycle number 999, 1999, 2999 ... The strange
numbers are because of the cycles being numbered from 0.

For instance, to extract the second configuration saved, at cycle 1999 in your
above simulation, you would write: \n

\verbatim
prompt> $profasi_dir/bin/extract_snapshot -c 2999 -o 2999.pdb n0/traj
\endverbatim

Notice that if your run took 5 segments to finish, there will be 5 conf.dataXXXX
files in your n0 directory. Read what is in the "traj" file. Verify that no
matter which conf.data... file contains the given cycle according to the traj
file, extract_snapshot finds it for you. The traj file presents an unbroken view
of the multi-segment run.

\subsection rtfl Run-time history file "rt"
The simulation writes a snapshot of different system properties to a "run-time
 history" file, almost always named "rt" in ProFASi applications. This file
 contains many columns: MC-time, temperature, total energy, different energy
 terms ... The MC time is measured in Monte Carlo cycles. The option "-nrt" sets
 the interval in MC time between the snapshots written to the rt file. The output
 directory n0, also contains a file called "rtkey", with information on the
 different columns in the rt file.

\subsection hisfls Histogram output files
The output directory "n0" also contains a number of files with names his_Etot,
his_Bias, his_HBMM etc. These are the histograms of different measurements
stored in the rt file. The histograms are filled after every sweep whereas the
rt file contains a snapshot every 1000 sweeps or so. So, the histograms have
more data. The format for storage in the histograms is obvious: some necessary
comment lines beginning with the character '#' followed by data points in x,y
pairs. (But this changes a bit when there are many temperatures! You will see
examples in tutorials on simulated or parallel tempering.)

Try plotting these histograms in gnuplot. You will see that all the
histograms show the data points well covered, and that they do not waste bins to
the left and right. This will always be the case with ProFASi version 1.1, and
later as the ProFASi histogram class AdaptiveHis adjusts its range to suit the
data, without losing any statistics. So, the user can get away without setting any
range.

The files with extensions ".profile", like "HelixContent.profile", store the
residue wise properties as a function of temperature. In this example, there is
one temperature. So, the data will just be one row with a number for each
residue for which such a property can be calculated. The HelixContent, for
instance, will skip the two terminal residues for which the Ramachandran angles
are not both defined.

\subsection avgfl The averages file
The output directory also contains a file called "averages", containing the
temperature wise averages of all measurements made during the run. For this
example, with a canonical Monte Carlo simulation, there is only one temperature.
The average file contains the index of the temperature, "0" in this case, and
the average values. For other simulation methods, like simulated annealing,
simulated tempering and parallel tempering, the averages file will contain
entries for temperature indices 0,1,2 ...

\subsection logfl The ProFASi logfile and other output files
There are a few other files in the output directory which are not related to the
measurements. The "logfile" contains log messages from different ProFASi classes
printed during the run. Browsing the log file often gives useful information
about what was happening during the run. It is especially useful for debugging.
The amount of messages printed in the log file can be controlled by setting the
log level for a run. If you run the above example with option "-ll 200", you
will see an explosion of new messages in the log file. Run it with "-ll 3" and
very little log messages will be printed.

That leaves the file, "updates.stats", which contains information about how often
 individual conformational update types of ProFASi were used, and how often they
 were accepted.

\section usingmimiqa Finding RMSDs to the native state
As mentioned, the native structure of the peptide used for this tutorial consists
of a simple alpha helix. How similar (or dissimilar) are the minimum
energy and other structures you extracted during the tutorial to the
experimentally determined structure ? There is a tool in ProFASi, called MimiqA,
to compare two protein structures in PDB files with the similarity measure
called, (minimized) Root Mean Square Deviation, RMSD. The peptide sequence used
above has a PDB id 1WFA. If you don't have the PDB file for the peptide, you can
download it from http://www.rcsb.org/pdb/home/home.do . Save the structure in the
run directory for this tutorial.

Now, to calculate the backbone RMSD, you should use the following command:
\verbatim
prompt> $profasi_dir/bin/mimiqa --using "+BB" 1WFA.pdb n0/minen.pdb
\endverbatim
For C-alpha RMSD, \n
\verbatim
prompt> $profasi_dir/bin/mimiqa 1WFA.pdb n0/minen.pdb -u "+CA"
\endverbatim
For RMSD over backbone and C-beta atoms, and restricting the calculation to
residues 2 -- 19 of the model 1 in the PDB file,
\verbatim
prompt> $profasi_dir/bin/mimiqa 1WFA.pdb:1:A,2,19 n0/minen.pdb --using "+BB+CB"
\endverbatim

The last example illustrates what happens when you use truncated residue ranges.
The simulated sequence has 37 residues while the selection 1WFA.pdb:1:A,2,19 has
18 residues. So, the sequences are not identical. ProFASi then performs a
sequence alignment before finding the RMSD. If a PDB file has less models than
the model id given on command line, the first model is used. To see what residues
got aligned with what, take a look at the log file that mimiqa generates:\n
\c prompt> less .profasi_logfile \n\n
This log file also contains an atom by atom comparison list. Convince yourself
that the order of atoms inside a residue in the two PDB files being compared
does not matter! For instance, the ring carbon atoms in Phenyl-alanine or
Tyrosine or Tryptophan residues appear in a different order in ProFASi than in
most files in the PDB. But you will see that CE1 will be compared to CE1, CE2
with CE2 etc. Of course, to see this you have to run mimiqa with the filter
"-u +HV" (HV for heavy atoms) in place of "-u +BB" etc.

MimiqA provides a great flexibility in calculating RMSD on the command line, and
can be used quite independently of the rest of ProFASi, in connection with data
analysis. The documentation of MimiqA has more details: \ref mimiqa

We remind the user that, here we are only introducing the tools of ProFASi. For
such simple helical peptides, it is not uncommon with ProFASi to end up with a
structure with a very low RMSD within 30 seconds of starting from a random state.
But neither such a thing in itself constitutes folding success nor the failure of
this example simulation to get a low RMSD constitutes a failure for the force
field. For a fair assessment, one should use runs lasting much longer and
typically with more sophisticated Monte Carlo methods, like simulated or parallel
tempering as will be explained in later tutorials.

\section animations What about an animation of the simulation in this tutorial ?
We end this tutorial with an explanation of how to make animations of a
particular part of the folding trajectory made with ProFASi. The utility we use
for this purpose is called extract_props . It can run through a specified range
of stored configurations and perform some operation on them. For instance,

\verbatim
prompt> $profasi_dir/bin/extract_props -op tmp_ -rt --start 999 --end 9999 n0/traj
\endverbatim

will recreate the run-time history of your run from the binary configuration
file, and create a new rt file called "tmp_rt". It should contain the same
numbers as the original rt file. Although it might appear to be a useless thing
to do at the moment, as you continue to use the program, sooner or later you will
come to the situation where you want to regenerate a run-time history based on
stored configurations.

Anyway, the program extract_props can also browse the configurations and create a
bunch of PDB files corresponding to those configurations. \n

\verbatim
prompt> $profasi_dir/bin/extract_props -op tmp_ -pdb -a 999 -z 9999 n0/traj
\endverbatim
This should create a series of PDB files, tmp_frame_1.pdb, tmp_frame_2.pdb etc.,
corresponding to the state of the molecule at the cycles where configurations
were saved. You can now proceed to render these snapshots with a PDB viewer as
image files and then combining the images into a movie with a program like
<a href="http://en.wikipedia.org/wiki/Mencoder">mencoder</a>, which is a
part of the <a href="http://www.mplayerhq.hu/design7/news.html">MPlayer</a>
package.

But you can also generate a single PDB file with the different configurations
along the trajectory being saved as different models in the file: \n

\verbatim
prompt> $profasi_dir/bin/extract_props -op tmp_ -pdb --single_file n0/traj
\endverbatim

This will create a single file, tmp_frames.pdb, with the snapshots at 999, 1999
... 9999 rendered as models 1, 2 ... 10. Omitting the "--start" and "--end"
options as above will result in the use of the entire range of available cycles.
You can load this PDB file in PyMol, and just "play" it! The PDB viewer will
display one structure after the other quickly, like a coarse grained animation.

The animation would be very "unsmooth" with large changes in conformation between
the frames. This is to be expected, as the configurations are 1000 Monte Carlo
cycles appart, which means there could be some 139000 Monte Carlo updates between
them, and we know that even a single Monte Carlo update can make large changes to
the conformation of a protein chain. But you can generate a slightly more refined
version of the animation, by re-running your simulation with -iconf 10, so that
you end up with 1000 snapshots. You can force the program to produce exactly the
same trajectory by using the random number seed option "-rs" as described above.
Only this time, it will save configurations more often. Then you can create a new
1000 model PDB file containing snapshots of the system every 10 Monte Carlo
cycles, instead of 1000.

\section nextch What next ?
This tutorial was only to get familiar with ProFASi. You will most probably not
run a canonical Monte Carlo simulation for a few minutes for a research project.
But now you can proceed to the next important component of ProFASi simulation
programs: the settings file. Using a settings file offers a somewhat more
convenient way of fine-tuning the behaviour of the simulation. That is the subject
of \ref tutorial_2.

*/