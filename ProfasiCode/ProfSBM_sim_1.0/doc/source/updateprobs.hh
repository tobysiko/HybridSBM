#ifndef UPDATEPROBS_HH
#define UPDATEPROBS_HH
/**
\page update_probs Specifying relative probability of conformation updates

As of now, PROFASI provides two  ways of fine tuning the conformational
updates. The first is a probability for choosing a conformational update, and
that will be described in this page. The other is the probability of using a
certain degree of freedom, after an update has been selected. That procedure
is described in \ref config_update. Assignment of relative probabilities as
described here, could also be done (most probably more elegantly) using the XML
configuration file, but that part has not been programmed yet.
<div style="color:red"> Consider this as historical baggage of the program which
will be removed in the near future, in favour of the much more elegant
mechanisms afforded by the PROFASI XML modules.</div>

At present there are 5 conformation updates implemented in PROFASI : Rot, Pivot,
BGS, Rotation and Translation. Rot acts on side chain DOFs alone. Pivot and BGS
change backbone DOFs. Pivot makes a random change of one selected angle. BGS
takes 4 consecutive residues and makes a concerted rotation to perform an
\e almost local deformation of the chain. Rotation and Translation are rigid
body updates operating on whole chains. During a Monte Carlo move, PROFASI
simulations using the MC class or its derivatives first select one of the 5
updates and perform them. Then they go on to accept or reject the move with
a Metropolis like criterion. The relative probabilities of the 5 updates
can be controlled by putting a file called "updateprob.dat" in the same
directory as the settings file. The contents of this file will be described
below, but first note that it is also possible to just ignore this file
all-together. If one starts simulations with no updateprob.dat file in the run
directory, the program analyses the system under study and automatically assigns
quite reasonable probabilities for the updates at different temperatures. The
automatic probabilities are based on proportions of backbone degrees of freedom,
side-chain DOF and rigid-body DOF in the system. Most of the time, you could use
the automatic probabilities, and worry about other things.

Now, if you decide on the probabilities using other criteria and want PROFASI to
use them, you can do so using the "updateprob.dat" file. Here is the syntax:

<ul>
<li> The file should begin with 3 lines with keywords "#Updates", "#T_indices",
and "#Columns". The rest of the lines are the actual probabilities ordered in a
matrix form.</li>
<li> The syntax and usage of the key words are:
<ol>
       <li> #Columns This specifies the orientation of your matrix, whether to
interpret different columns as information about different temperatures or as
different updates. It can have only two reasonable forms: <i>#Columns
Temperatures</i> or <i>#Columns Updates</i></li>
       <li> #Updates This is used to specify the updates for which the file
contains probability information. For example:\n
#Updates BGS Pivot Rotation Translation Rot

This tells the program that the file contains probabilities for the 5 updates in
that order.</li>
       <li> #T_indices This specifies the order in which the temperature indices
are to be interpreted for this file. Normally PROFASI orders the temperatures in
a decreasing order, so that a temperature index of 0 corresponds to the highest
temperature. But you might find it more convenient to think in terms of
increasing temperatures and write your probability matrix, where increasing row
(or column) index corresponds to increasing temperature. You would then ask the
program to interpret the temperatures in reverse. For instance:\n
#T_indices 7 - 0\n
PROFASI will then interpret your first row or column as temperature 7 and so on.
</li>
</ol>
</li>
<li> All other lines in the file are interpreted as a space separated list of
numbers.</li>
<li> The number of fields in the lines need not be constant. If all lines have
5 numbers and one line has 3, then the last of those 3 numbers will be copied 2
more times so that all lines will have the same number of entries.</li>
<li> If an update is specified which is not known by the program or is not
currently being used, the row or column corresponding to that update will not be
used.</li>
<li> If the "#T_indices" command is used to specify order for N temperatures but
the program uses ntmp < N temperatures, only ntmp values corresponding to the
highest temperatures in your file will be read, i.e. indices 0 to ntmp.If on the
other hand, ntmp > N, the probabilities for the unspecified temperatures will be
assumed to be the same as for the last specified temperature.</li>
<li> The probability values for any temperature are normalized in the program.
So, as long as at least one update receives a non-zero value for a given
temperature, the specified probability values do not need to be normalized in
the updateprobs.dat file.</li>

</ul>

<b>Example 1</b>\n
First a case where the file "updateprobs.dat" contains the following lines
\verbatim
#Columns Updates
#Updates Rot BGS Pivot Unknown Translation Rotation
#T_indices 7 - 0
0.7    0.25    0.025    0.4    0.0    0.025
0.7    0.15    0.075    0.4    0.0    0.075
0.7    0.05    0.125    0.4    0.0    0.125
\endverbatim
The following probabilities will be assgned:\n
Rot:(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)\n
Pivot:(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.075, 0.025)\n
BGS:(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.15, 0.25)\n
Translation:(0, 0, 0, 0, 0, 0, 0, 0)\n
Rotation:(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.075, 0.025)\n

Notice that only 3 columns corresponding to 3 temperatures are given, and since
the temperature indices are intreptreted 7 - 0, it is inferred that temperatures
7,6 and 5 have been specified and the remaining not. So, the values for
temperatures 0 though 4 are copies of the values read for temperature 5. Notice
also that the column for the update called Unknown here is completely skipped,
as the update is not known.


<b>Example 2</b>\n
\verbatim
#Columns Temperatures
#Updates BGS Pivot Rot
#T_indices 0 - 7
0.05    0.05    0.1    0.15    0.2    0.25
0.25    0.25    0.2    0.15    0.1    0.05
0.7
\endverbatim
The following probabilities will be assgned:\n
Rot:(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)\n
Pivot:(0.25, 0.25, 0.2, 0.15, 0.1, 0.05, 0.05, 0.05)\n
BGS:(0.05, 0.05, 0.1, 0.15, 0.2, 0.25, 0.25, 0.25)\n

The first row of the probability matrix is used for BGS, second for Pivot and
third for Rot. The Rot line has equal entries for all temperatures, so it is
sufficient to specify it just once.
*/
#endif // UPDATEPROBS_HH
