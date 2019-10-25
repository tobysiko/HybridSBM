/**
* \page tutorial_2 Tutorial 2: Using the PROFASI settings file
In tutorial 1, we used command line options to control the behaviour of various
PROFASI programs. It is not convenient to pass 25 command line options to fine
tune a complex parallel tempering run. Therefore, PROFASI simulation programs,
such as BasicMCRun, SimTempRun and ParTempRun, also read commands from a
"settings" file in the run directory, in addition to processing the command line.
The settings file has a fixed name: "settings.cnf", and is located in the run
directory. Most of the options you passed to the BasicMCRun program in tutorial 1,
can be set through the settings file, but it is more convenient and tidy to work
with a settings file.

\section example_run An example simulation: Requesting run time measurements
To explain the use of the settings file, we will rerun the simulation of tutorial 1.
Only this time, we will add two new columns to the run-time history file "rt" to
record the radius of gyration and the RMSD to the native state. We will also
obtain histograms of these two properties for the structures seen during the run.
Since we will calculate RMSD during the run, the file containing the native state,
1WFA.pdb, is needed. So, create a new run directory for this tutorial and copy
1WFA.pdb into it. Now using a text editor create a new file called "settings.cnf"
in the run directory, and type in the following text into it...\n\n


\verbatim
new_obs Rg rg
new_obs ProteinRMSD bbcb using +BB+CB; struc1 1WFA.pdb:1:A,2,19 ; struc2 $::A
\endverbatim


Save the file, and run the simulation exactly like in tutorial 1:
\verbatim
prompt> $profasi_dir/bin/BasicMCRun \
-ac 1 "* DTASDAAAAAALTAANAKAAAELTAANAAAAAAATAR * NH2" \
-T "300 Kelvin" -ncyc 10000 -nrt 1000 -iavg 1000
\endverbatim

The two commands in the settings file created two new "Observables" for the run.
The first argument to the "new_obs" command is the name of a PROFASI Observable
class that implements a particular kind of measurement. The second is an alias
given by the user to the newly created observable. Notice that different kinds
of measurements require different set up information from the user. RMSD needs
to know what structure to compare with, what kind of atoms to use etc. Radius of
gyration does not need a comparison structure, although even that might be
restricted to a range of residues. It would be rather cumbersome to do all this
on the command line for a simulation program. Thus the settings file.

You should now check the "rtkey" file in the output directory and see that there
are indeed two more observables listed by the aliases "rg" and "bbcb".
Correspondingly, the rt file will contain additional columns compared to
tutorial 1. Check also that there are two more histogram files for these two
observables, and that the program chose reasonable ranges for them.

To choose and fine tune measurements to be performed during the run is the main
purpose of the settings file. To learn about how to set up other observables
using the settings file, you should consult \ref settings_obs .

\section set_args Relation between settings file commands and the command line
In PROFASI 1.5, the long versions of the command line options also work as
settings file instructions, so that there is only one kind of instructions to
learn. For instance, the \e -T "300 Kelvin" option above is the short form of
\e --temperature "300 Kelvin" . \e -ac is the short form of \e --add_chain . One
could use either the long or the short forms of the commands on the command line.
The documentation of the commands will list both forms. The long version without
the initial "--" is the equivalent settings file instruction. The only difference
is that multiple word arguments for instructions, like "300 Kelvin", need to be
in quotes on the command line, and <b> must be </b> without quotes in the
settings file. The add_chain command has two arguments: the number of chains and
the sequence. It looks for two arguments. That's why the 1 is not inside the
quotes containing the sequencein the above example. In the settings file, since
each instruction starts in a new line, the quotes are not needed. The first
example in the next section will use these two in addition to a lot of other
instructions in the settings file.

\note Note that in PROFASI 1.5, if a settings file instruction does not
fit in one line and needs to be continued across many lines, the continuation
character "\" needs to be used at the end of the incomplete lines. This is
also demonstrated in the example below, although it is not strictly necessary
there.

\section settings_all Settings file without any command line options
Since we have to use a settings file for a part of the task, it would be
convenient if there was a possibility to set up the run entirely using it,
without the command line options for sequence etc. This is possible, and this
way of using the program is more tidy. Consider the following settings file for
a run: \n

\anchor tut2ex2
<b>Example 2</b>
\verbatim
log_level 10
add_chain 2 <*KFFE AAAK KFFE*>
add_chain 2 <*KFFE YNGK \
              KFFE*>
box_length 60
temperature 290 Kelvin
rt_write_freq 1000
conf_write_freq 1000
steps_per_cycle 200
num_cycles 10000
new_obs Rg rg0 of_chain 0
new_obs Rg rg1 of_chain 1
new_obs Rg rg2 of_chain 2
new_obs Rg rg3 of_chain 3
new_obs ProteinRMSD rmsd using +HV ; struc1 $::A ; struc2 $::B
\endverbatim

If you save the above in a file called "settings.cnf" and just type
\verbatim
$profasi_dir/bin/BasicMCRun
\endverbatim

you will start a canonical Monte Carlo simulation with 4 peptide chains in a
60 &Aring; periodic box. Of the 4 peptide chains, 2 have one sequence and the
rest have another.

A lot of these commands in the settings file are simply the long versions of
familiar command line options. The "lcyc" option used to set the length of
the Monte Carlo cycle is now spelt out as "steps_per_cycle". Here we explicitly
set the size of cycle to 200 elementary MC steps. By default in PROFASI 1.5, a
cycle consists of as many steps as the number of degrees of freedom.

We have also added 4 different radius of gyration observables measuring the
property of the individual chains. We also introduce a novelty: RMSD between
two "live" structures in the simulation. The observable rmsd will measure the
root mean square difference between the first two chains added during the course
of the simulation. It's a measure that will tell us if the two chains adopt
similar structures at the same time.

What happens when an option given in the settings file contradicts one given as
a command line argument ? For instance, using the above settings file, you can
start the program with
\verbatim
$profasi_dir/bin/BasicMCRun -ncyc 20000
\endverbatim
In such a case, the value given on the command line has priority. So, in this
case, the program will run for 20000 Monte Carlo cycles and not 10000. If you
have set up a long run, but want to check that your set up makes sense, you can
override the number of cycles with a small value for a brief check run without
changing the settings file. But note that if you run
\verbatim
$profasi_dir/bin/BasicMCRun -ac 10 "*KLVFFAE*" -ncyc 20000
\endverbatim
you do not replace the simulated system with 10 KLVFFAE peptides. You create a
system with 10 KLVFFAE, 2 KFFEAAAKKFFE and 2 KFFEYNGKKFFE peptides instead. It
makes perfect sense to have many chains of different kinds in the system, and
therefore, having sequence information on the command line and the settings file
are not really "contradictory".

More details on how to use the settings file can be found in
\ref settings_population, \ref settings_random_number, \ref settings_mc,
\ref gmc_opts and \ref settings_obs.

In the next few tutorials we will introduce real MC simulations with PROFASI.

\li \ref tutorial_1
\li \ref tutorial_3
\li \ref toc

*/

