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

#include "Contact.hh"

namespace prf
{
    Contact::Contact(): first(0),second(0) {}

    Contact::~Contact() {}

    Contact::Contact(int i1,int i2) {first=i1;second=i2;}

    Contact::Contact(const Contact &c) {first=c.first;second=c.second;}

    Contact &Contact::operator=(const Contact &c)
    {
        if (this!=&c) {first=c.first;second=c.second;}

        return *this;
    }
}

