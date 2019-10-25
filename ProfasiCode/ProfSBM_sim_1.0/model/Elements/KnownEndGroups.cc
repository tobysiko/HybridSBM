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

#include "KnownEndGroups.hh"

namespace prf
{
    EndGroup * new_EG_object(prf::OneLetterCode cd)
    {
        switch (cd) {
            case ACE:
                return new Acetyl();
            case SUC:
                return new Succinyl();
            case NH2:
                return new Amide();
            case NME:
                return new NMethyl();
            case VOIDEG:
                return new VoidEG();
            default:
                return NULL;
        };
    }
}
