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

#ifndef SeqBuild_HH
#define SeqBuild_HH
#include <string>
#include "../Elements/GroupLib.hh"
#include <vector>

namespace prf
{

    class SeqBuild
    {
    public:
        SeqBuild();
        ~SeqBuild();
        void parse(std::string inpstr,std::vector<OneLetterCode> &sq);
        std::string make_string(const std::vector<OneLetterCode> &sq,
                                char begn='<', char nd='>');
        std::string next_word(int &ist);
        char next_char(int &ist);
        OneLetterCode next_code(int &ist);
        inline void word_mode() {charmode=false;}

        inline void letter_mode() {charmode=true;}

    private:
        bool charmode;
        std::string workstr;
    };
}

/**
\page seq_input Accepted sequence input format in PROFASI

The sequence information can be given in spelled out residue names (if it is a
single word), 3 letter codes or one letter codes. By default, sequence reading starts
in the word mode: either 3 letter codes or full names. If the character
"*" (asterix) is encountered, the reading mode toggles between word and one letter code
modes. So, <i>*VAL VAL*</i> is the sequence V-A-L-V-A-L, and
<i>VAL VAL</i> is the sequence V-V. If there are capping groups, like Acetyl (ACE),
they can be specified with their names or 3-letter codes.

Other possible capping groups are C-terminal amide (NH2) or N-methyl (NME) and the
N-terminal Succinylic acid (SUC) and the "VoidEG". The VoidEG (written like that
in the sequence for PROFASI) takes the place of a capping group, but does not
add anything there. This forces the program to start the chain at the N-terminus
with a NH-Calpha... instead of a charged N terminus with 3 hydrogens. In other
words, the residue at the N-terminal is treated as if it were in some other
location on the chain. VoidEG can be attached also at the C-terminus for a
similar effect.

As a contrived example, let's consider a peptide
Acetyl-ACDEFGHI-(D-proline)-KLMNPQRSTVWY. Let's also imagine that for some reason
we want the C-terminal Y to behave as if it was in the middle of a chain: i.e.,
no COO(-) ending. In the settings file for simulation programs, we would enter
this sequence as:

\verbatim
add_chain 1 < ACE *ACDEFGHI *DPR* KLMNPQRSTVWY *VoidEG>
\endverbatim

We start writing in the "word mode" and enter ACE for Acetyl. Then with the "*"
we toggle to the "letter mode" and enter the sequence Alanine (A), Cysteine (C),
Aspartic acid (D) ... up to Isoleucine (I). Then we need D-proline which does
not have a single character code like natural amino acids. So, we use the "*" to
toggle to the word mode and enter its 3-letter code DPR, and toggle back to the
letter mode for the regular amino acids. Then we continue with Lysine (K),
Leucine (L) ... until Tyrosine (Y). We switch to the word mode again and name
the special capping group VoidEG which fakes the presence of a capping group
for the chain constructor routines without adding any atom.

*/

#endif
