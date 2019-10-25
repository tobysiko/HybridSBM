//Contents of this file are required only for documentation

/**
 \page installation Installation instructions for ProFASi version 1.5
You will need a good fairly recent C++ compiler, and the GNU "make" utility.

PROFASI comes with a Makefile to manage the compilation. If your compiler is one
of GNU g++, Intel icc, IBM xlC or Portland Group's pgCC compilers, you can
start compilation directly by running "make" with an argument "CC=your_compiler"
in the main PROFASI directory. For instance, if your compiler is called "icc",
you would run, <br>

\verbatim
make clean
make CC=icc
make CC=icc parallel
make install
\endverbatim

The install step here is only a "local install". So, have no fear and do it! It
only creates an "include" and a "bin" directory under PROFASI/installs, and
installs the program executables and headers in them. You do not need root
priviledges for that. If you are using git as a version control system, make
install creates a sub-directory with the name of the current git branch, and
installs into it. For instance, if you are in the "master" branch, and you
compile and install, you will have a directory PROFASI/installs/master/bin with
all executables and PROFASI/installs/master/include with all headers.

The MPI dependent parts of the program are only compiled in the make ___ parallel
step. It is possible to skip that step, and use all the serial parts of the
program, including all the analysis tools and even most of the simulation programs
such as BasicMCRun, SimAnnealRun and SimTempRun. There are parallel versions of
these programs which will not be built if you skip the make __ parallel step. If
you have MPICH set up for a certain compiler in your system, you should pass that
compiler to the CC flag.

Some good optimisation options have been chosen for each of these compilers in
the Makefile. If your compiler is not one of these, the optimisation and other
options it takes might be different, and you will need to change the CC, CFLAGS
and possibly MPICC and MPICFLAGS options in the file \em Makefile.compilers
before you can compile.

On some systems, you might need to invoke the GNU make utility with the command
"gmake".

The make process enters the subdirectories <i>model</i> and <i>app</i> and runs
make inside them. If compilation succeeds, you should have a library
libprofasi_(VERSION).a created in the directory <i>model/lib</i>, a set of
executable files in the directory <i>app/bin</i>.

If the library fails to build, or if there is a problem with linking so
that the executables cannot be built, and you can pinpoint the problem, we
would like to hear about it.

The present release has been tested with the following compilers: <br>
<ul>
<li>gcc version 4.4.3 </li>
<li>Intel(R) C++ compiler version 11.1</li>
<li>Portlan Group's compiler suite version 10.1</li>
<li>IBM(R) XL C/C++ Enterprise Edition V8.0</li>
</ul>
It is not possible for us to guarantee that the code will compile with older
compilers. The compilers tend to get better with time. So, if PROFASI compiles
with the latest compilers, we are satisfied. Since there is an open source C++
compiler with a  reasonably high standard (gcc), we do not offer any help in
compiling our code on ancient compilers.  If, on the other hand, you notice that
PROFASI fails to compile with a newly released version of a reputable compiler,
please let us know.

We have no possibility for ensuring compilation and execution in a Microsoft
Windows environment.

*/

