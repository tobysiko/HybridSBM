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

#include "NaturalAminoAcids.hh"

namespace prf
{
    AminoAcid * new_AA_object(prf::OneLetterCode cd)
    {
        switch (cd) {
            case G:
                return new Glycine();
            case A:
                return new Alanine();
            case V:
                return new Valine();
            case L:
                return new Leucine();
            case I:
                return new Isoleucine();
            case S:
                return new Serine();
            case T:
                return new Threonine();
            case C:
                return new Cysteine();
            case M:
                return new Methionine();
            case P:
                return new Proline();
            case DPR:
                return new Proline(DEX);
            case D:
                return new Aspartic_Acid();
            case N:
                return new Asparagine();
            case E:
                return new Glutamic_Acid();
            case Q:
                return new Glutamine();
            case K:
                return new Lysine();
            case R:
                return new Arginine();
            case H:
                return new Histidine();
            case F:
                return new Phenylalanine();
            case Y:
                return new Tyrosine();
            case W:
                return new Tryptophan();
            default:
                return NULL;
        };
    }
}
