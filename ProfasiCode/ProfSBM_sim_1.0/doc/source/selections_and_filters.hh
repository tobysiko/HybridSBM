/**
\page prf_sel_fils ProFASi PDB handling: Selections and Filters
\section selections Selecting parts of one or more chains
A selection is a set of residues upon which some operation is to be performed. 
In its most simple form, it is a contiguous range of residues, say, the first
to the tenth residue of a chain( with a chain label "A"). This selection can
be conveniently written as "A,1,10". A simple selection in PROFASI is always 
of this form: a chain label character and two residue indices, separated 
by commas. 

A selection may consist of non-contiguous regions. Say, the 
residues 1-10 like above, and residues 26-35. Such compound selections are
written by joining simple selections with a colon (":") mark in between. 
In this case, it would be "A,1,10:A,26,35". A compound selection can have 
segments of different chains together. For instance, "A,1,10:B,1,10" means the
object formed by taking the residues 1-10 of chains A and B together. A selection
of an entire chain is done by omiting the residue range. That is, a selection
string "A", means the entire chain A. Also, "A,10" means residues 10 through 
last, and "A,,19" means residues first available through to 19. 

In residue selections, the range indicated is assumed to be inclusive on both
ends. That is, "A,1,5" means residues 1,2,3,4,5 in chain A. 

The residue indices in the range are strings rather than integers. So, if residues
10,11, and 12 are missing in a PDB file, a selection "A,8-15" selects residues 8,9,
13,14,15. If there are two extra residues between 10 and 11, called 10A and 10B
in the PDB file, the same selection will pick 8,9,10,10A,10B,11,12,13,14,15.  

\section filters Filters: choosing a limited set of atoms for an operation
It is often necessary to select a certain kind of atom, like all C_alpha atoms, 
and do something with them. The concept of "selections", as described in 
\ref selections , only applies to ranges of residues. The selection of atoms
based on their attributes is done with "filters". 

A filter is a criterion that decides whether or not an atom is to be included
in the selection. 

In an application program, the interface to the filters is intuitive. There are 
inclusion filters and exclusion filters. The simplest kind of filters work with 
one single attribute of the atom. It is best to explain using a few examples...

\li "+@C" will select all carbon atoms. "-@S" will exclude sulfur atoms.

\li "+CA" will select all atoms labeled "CA", i.e., the C_alpha atoms. Similarly
"+CB", "+N", "+OG" ... any valid label for an atom in a PDB file.

\li "-%PRO" will exclude any atom that is in a Proline residue.

\li "+BB" is a shorthand for a filter to accept an atom if it is a backbone
N, CA or C. "+HV" is an alias for "-@H", that is, a filter to accept any atom
that is not hydrogen. 

Filters can be combined. To select all C_alpha and C_beta atoms, one puts
the inclusion filters together: "+CA+CB". To select all heavy atoms but not
sulfur atoms, one uses "+HV-@S". The facility to combine any number of
inclusion and exclusion filters  provides a way to make very complicated 
selection criteria on the fly. 

In PROFASI, filters are implemented as "predicates" of prf::AtomDescriptor 
objects. A predicate is a functional that acts on one type of object and returns 
either true or false. PROFASI contains a composite capable predicate class in
\ref CompositePredicate . That is, if you have a predicate, p, to select C_alpha 
atoms; and another, q,  to select C_beta atoms, you can combine them with p||q, 
and create a new predicate to select an atom if it is either C_alpha or C_beta. 
Predicates can be combined with "or", "and" and "xor" operations, and they can be 
negated. This allows very complex inclusion or exclusion criteria to be built up. 

\sa \ref mimiqa
*/
