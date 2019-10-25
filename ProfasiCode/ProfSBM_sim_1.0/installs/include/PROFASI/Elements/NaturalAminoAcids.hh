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

#ifndef NATURAL_AMINO_ACIDS_HH
#define NATURAL_AMINO_ACIDS_HH

#include "../AA/glycine.hh"
#include "../AA/alanine.hh"
#include "../AA/valine.hh"
#include "../AA/leucine.hh"
#include "../AA/isoleucine.hh"
#include "../AA/serine.hh"
#include "../AA/threonine.hh"
#include "../AA/methionine.hh"
#include "../AA/cysteine.hh"
#include "../AA/proline.hh"
#include "../AA/aspartic_acid.hh"
#include "../AA/asparagine.hh"
#include "../AA/glutamic_acid.hh"
#include "../AA/glutamine.hh"
#include "../AA/lysine.hh"
#include "../AA/arginine.hh"
#include "../AA/histidine.hh"
#include "../AA/phenylalanine.hh"
#include "../AA/tyrosine.hh"
#include "../AA/tryptophan.hh"

namespace prf
{
    AminoAcid * new_AA_object(prf::OneLetterCode cd);
}

#endif
