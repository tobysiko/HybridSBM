/**
\page tutorial_3 Tutorial 3: Simulated annealing
In simulated annealing, the temperature for the canonical Monte Carlo run is
gradually lowered, starting from some high temperature. We illustrate simulated
annealing in PROFASI and its peculiarities by taking an \f$\alpha\f$-helical
mini-protein: the Tryptophan cage 1L2Y. Instead of BasicMCRun, we now run
SimAnnealRun: \n\n

\verbatim
prompt> $profasi_dir/bin/SimAnnealRun -ac 1 "*NLYIQWLKDGGPSSGRPPPS*" \
-ncyc 20000 -nrt 500 -tmax "374 Kelvin" -tmin "274 Kelvin" -ntmp 10 -ncT 1000
\endverbatim

First, notice that the usage is very similar to that of BasicMCRun. This is the
case with all simulation programs in PROFASI. Technically, SimAnnealRun,
SimTempRun, ParTempRun and WLRun inherit from BasicMCRun. So, their interfaces
will seem similar. Some options for BasicMCRun, like the <i>--add_chains</i>
option, make perfect sense here, and are simply inheritted from BasicMCRun. Some
options, like <i>--temperature</i> make no sense here and are therefore ignored.
Notice that in the above, we use options \e -tmax and \e -tmin to set a range
of temperatures instead of specifying one temperature. \e -ntmp means the number
of temperatures going from the maximum to the minimum. The option \e ncT specifies
how many Monte Carlo cycles are to be performed at each temperature. The
temperatures are lowered by a constant factor every ncT Monte Carlo sweeps so
that we go from 374 K to 274 K in 9 changes. You can of course choose to have
ncT=1 and 100000 temperatures in your simulation, although this is not
recommended as the averages will be printed for all temperatures in the
simulation, leading to some huge output files.

\section tutorial_3_output Output files
Now, let's examine the output of the program. On the standard output, the
simulated annealing program writes down what temperature values are used, how
many Monte Carlo sweeps are performed at each temperature, and the probability
for various conformational updates to be used at each temperature. Now, open the
run-time history file, rt. You will notice that the second column, which always
used to be 0 for the BasicMCRun examples, now has non-zero values. The value in
this column is the "temperature index", and indicates which temperature the
system is at, at a certain MC time. The actual temperature can be looked up in
the standard output list of temperatures, as mentioned above, or in a special
file called <i>temperature.info</i> in the output directory. If you plot the
column 2 of rt vs column 1, you should see a saw-tooth pattern. The index
increases linearly from 0 to it's maximum value (the number of temperatures-1),
and then goes to 0. Since PROFASI orders temperatures by
\f$\beta = \frac{1}{k_B T} \f$, increasing temperature indices normally
correspond to decreasing temperatures.

The "averages" file, which contained a single row of values for each measurement
in BasicMCRun, now has several rows for each measurement, one for each
temperature. The values are listed according to the temperature-index.

Histograms are turned off by default for simulated annealing in PROFASI although
you can enable them by run by doing both of the following:

\li use the option \e --histograms on the command line
\li put the number of thermalisation cycles to 0 with \e -ntherm \e 0, as the
program does not send any data to the histograms until the designated number of
thermalisation cycles are finished. Since simulated annealing spends a
pre-determined number of cycles at each temperature, it is better to set ntherm
to 0 and have the same number of points for each temperature.

If you do the above, you will notice that the generated histogram files look
different and have many columns. The first column is the value of the observable
in question, say x, and the subsequent columns are the dP/dx values at
temperature indexes 0, 1, 2, 3 ... Each histogram file contains data for all
temperatures. (Yet another thing to consider if you want to do simulated
annealing with 100000 temperatures!) The reader will learn more about PROFASI's
histograms in subsequent tutorials on simulated and parallel tempering.

Most other output files have similar meanings in SimAnnealRun as in BasicMCRun.

\section tutorial_3_finegrained Schema files: fine grained control over temperature decrease

It has often been suggested that it is a waste of computation time to have as
many MC cycles at the highest temperatures of simulated annealing as at the
lowest temperatures. This is not the place for us to argue in favour or against
such a claim. But we provide facilities in PROFASI to use virtually any scheme
for lowering temperatures in simulated annealing. This is done with a "schema"
file for simulated annealing. The following is an example schema file:
\verbatim
#temperature Kelvin
400 100
350 120
330 500
320 1000
315 3000
310 1000
300 1000
290 1000
280 1000
270 1000
\endverbatim
The first line tells the program that the temperature values are in Kelvin.
Alternatives would be "#temperature inverted" if you choose to input
\f$\beta = \frac{1}{k_B T} \f$ values, or you omit the "#temperature ___" line
altogether so that the temperature values will be assumed to be in internal
PROFASI units. The remaining lines are interpreted as ordered pairs with a
temperature and the number of sweeps at that temperature. So, using this schema,
the simulated annealing program will spend only 100 sweeps at the highest
temperature, 120 at 350 K, ... To use this scheme of lowering temperatures, save
the above example schema file (say, as abcd.sas), and invoke the SimAnnealRun as,

\verbatim
prompt> $profasi_dir/bin/SimAnnealRun -ac 1 "*NLYIQWLKDGGPSSGRPPPS*" -ncyc 20000 -nrt 500 --schema_file abcd.sas
\endverbatim

In case you use such a schema file, options tmin, tmax, ntmp and ncT are ignored.
All such information can be inferred from the schema file.

You can also run the simulated annealing program without these command line
options, but with a settings file like the following:

\verbatim
log_level 10
add_chain_pdb 1 1L2Y.pdb:1:A,1,20
box_length 125
tmin 200 Kelvin
tmax 400 Kelvin
ntmp 50
ncT 1000
nrt 1000
iconf 1000
ncycles 200000
new_obs Rg rg
new_obs ProteinRMSD bbcb using +BB+CB; struc1 1L2Y.pdb:1:A,2,19 ; struc2 $::A
\endverbatim

Note how we used a different way to input the sequence. Since we are calculating
RMSD, the PDB file is already in the run directory. So, why not extract the
sequence form there ?!

You can omit the commands tmin, tmax, ntmp and ncT and provide a schema file with
the command "schema_file" in the settings file just like we did on the command
line earlier.

\section tutorial_3_common_error A common error
We have often received error reports when the user tries something like this and
is puzzled why simulated annealing thinks it has nothing to do.
\verbatim
prompt> $profasi_dir/bin/SimAnnealRun -ac 1 "*NLYIQWLKDGGPSSGRPPPS*" \
-ncyc 20000 -nrt 500 -tmax "374 Kelvin" -tmin "274 Kelvin" -ntmp 50 -ncT 1000
\endverbatim

If you look more carefully, we are asking for 50 temperatures and 1000 MC
sweeps at each temperature. That requires 50000 sweeps, whereas we have given
the program only 20000 sweeps. Simulated annealing only performs an integral
number of annealing cycles from maximum to minimum temperature. It calculates
the number of annealing cycles it can complete and then adjusts the number of
MC cycles to match that. For instance, in the above example, if we had used
\e -ntmp \e 9 instead of 50, the standard output would have informed us that
the program was trying to get from cycle 0 to 18000 even if we asked for 20000.
In the case of the above example, that readjusted MC cycle count becomes 0, and
hence, nothing happens.

\section tutorial_3_parallel Running simulated annealing in parallel
You can also run several simulated annealing runs in parallel conveniently by
running the binary "SimAnnealRun.mex" (generated if you did a "make parallel"
when you compiled PROFASI).

\verbatim
prompt> mpirun -np 64 $profasi_dir/bin/SimAnnealRun.mex -ac 1 "*NLYIQWLKDGGPSSGRPPPS*" -ncyc 20000 -nrt 500 --schema_file abcd.sas
\endverbatim

There is no communication between different parallel ranks in simulated
annealing. SimAnnealRun.mex is simply a "farming" approach to exploit available
computing power in computing clusters to get more annealing cycles quicker.
Parallel tempering, covered in a later tutorial, is a genuinely parallel
simulation algorithm.

\section tutorial_3_closing Closing remarks
We remind the reader that simulated annealing is a global optimisation method.
Often it is useful to find the minimum energy conformation for a given system
in a force field, but it is not particularly useful to calculate the temperature
dependence of measurable quantities. Even for finding the global minimum of
energy in protein systems, if we start with no information about the ground
state, our experience indicates that replica exchange or parallel tempering works
noticeably better for all but the simplest of peptide systems.

If you use this method, you should make sure that there are a large number of
annealing cycles, going from the maximum to the minimum temperature, either by
having very long trajectories with many annealing cycles or a large number of
runs in parallel.

\li \ref tutorial_2
\li \ref tutorial_4
\li \ref toc
*/

