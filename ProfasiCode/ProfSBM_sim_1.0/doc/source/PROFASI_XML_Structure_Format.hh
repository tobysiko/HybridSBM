#ifndef PROFASI_XML_STRUCTURE_FORMAT_HH
#define PROFASI_XML_STRUCTURE_FORMAT_HH
/**
\page xmlstruc PROFASI's XML structure file format
While working with PROFASI, it is often convenient to use its own XML structure
format to store the structure of protein systems. There are many reasons for
introducing this:
\li The punch-card age conventions in the PDB format can sometimes be somewhat
invonvenient, for instance when one has more than 100 chains. The chain label
column with only one character can be a limitation, leading to chain labels
such as ',','?' or unprintable characters.
\li Although specifying 3 digits for the Cartessian coordinates is normally ok,
one can not save the precise state of the population in that way. For instance,
if the program was at a state S with energy E, and we save a PDB file. The energy
calculated for the PDB file will differ by a small but non-vanishing amount.
This is especially annoying if you are trying to develop a force field, and
can not be sure if the small differences you see are because of a change of your
implementation of the energy function or just because of the rounding off errors
in the PDB format.
\li Earlier versions of PROFASI have had a "textconf" format, in which only the
degrees of freedom are written to a file. This format was what the developers
used to restore the population to a stored structure precisely. But since
textconf does not contain sequence information, it can only be interpreted in
connection with a settings file. It is desirable to have something that stores
structures for PROFASI as accurately as textconf, but has the self-contained
character of a PDB file.

All the above issues are addressed by PROFASI's XML structure format. Here is
a sample:
\verbatim
<?xml version="1.0"?>
<structure>
<box_length>57.6999999999999957</box_length>
<energy>127.6096792213545825</energy>
<snapshot_time> 0 </snapshot_time>
<profasi_version>1.4.8</profasi_version>
<creation_time>
UTC 2011-Feb-01-17:44:59</creation_time>
<population>
<num_chains>1</num_chains>
<protein id="0">
<sequence>
  * RGKWTYNGITYEGR *
</sequence>
<global_coordinates>
-6.2660000000000000  -5.5970000000000004  -3.4409999999999998
-5.9421299751309578  -5.6627480501613841  -2.0188940261389439
-4.7606256705779826  -6.5864440712716199  -1.7714405175602406
</global_coordinates>
<group index="0" type="ARG">
<coordinates>
2.0942352533962136   2.2435194293934790   -3.1415926535897931
2.1414162734233773   -2.6370843727438729   1.5540134319459780
2.8388317673757530
</coordinates>
</group>
<group index="1" type="GLY">
<coordinates>
1.8356117390315918   0.8009091894546179   -3.1415926535897931

</coordinates>
</group>
<group index="2" type="LYS">
<coordinates>
-2.8358781065452923   2.3615227222406134   -3.1415926535897931
-2.3327310248966375   -2.7777380595288519   -1.8560468038439331
2.8623692260015430   -1.0468954991742636
</coordinates>
</group>

///// MANY MORE NODES /////

<group index="13" type="ARG">
<coordinates>
-1.7226977924681979   -0.2137330785775413   -3.1415926535897931
-2.3538889062619526   -2.7744724772572389   -2.1885416072649102
3.0158907181494810
</coordinates>
</group>
</protein>
</population>
</structure>
\endverbatim

The meaning of many of the fields does not need any explanation.
\li There is a  XML node tree for \e population, with \e protein nodes as
children.
\li Each protein has a \e sequence node containing its sequence, and a set of
global coordinates, storing the Cartessian coordinates of the first 3 backbone
atoms from the N-terminal.
\li The protein node has a bunch of residue nodes as children.
\li Each residue  node (called "group") has an "index" attribute, which
specifies where in the sequence the residue should be mapped. It also records
its 3 letter code in the "type" attribute.
\li The residue or group then has a "coordinates" child node containing a few
numbers which specify the conformation of that residue. The reason it is called
"coordinates" and not "dof" is that we include also the backbone omega angles
in this which are currently not degrees of freedom in PROFASI. The value of
the omega angle is read and used to decide if it is a CIS or a TRANS peptide
bond.
\li The coordinates node stores in this order: phi, psi, chi0, chi1 ... The
phi angle of Proline is stored, although it is not read in when the structure
is read.

The XML format can be displayed using a web browser like firefox and processed
with a plethora of different applications. White space is unimporant.

Since no visualisation program is likely to understand this format yet, we need
a way to convert the XML files to PDB files. This is how one does it:

\verbatim
prf_convert abc.xml abc.pdb
\endverbatim

The simulation programs can be asked to generate the minimum energy structures
in the XML format instad of the PDB format with the following command in the
settings file:

\verbatim
snapshot_format 2
\endverbatim

Otherwise, one can extract XML structure files from a concluded simulation using
extract_snapshot .

\sa \ref prf_convert, extract_snapshot

*/
#endif // PROFASI_XML_STRUCTURE_FORMAT_HH
