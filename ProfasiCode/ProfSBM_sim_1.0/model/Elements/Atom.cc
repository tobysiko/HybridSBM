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

#include "Atom.hh"

using namespace prf;

Atom::Atom() : my_type(hydrogen), pos(0) {}

Atom::Atom(int i) : my_type(hydrogen), pos(i) {}

Atom::Atom(int i, AtomKind tp) : my_type(tp), pos(i) {}

Atom::Atom(const Atom &atm) : my_type(atm.my_type), pos(atm.pos.Id()) {}

Atom::~Atom() {}

Atom & Atom::operator=(const Atom &ga)
{
    if (this != &ga) {
        my_type=ga.my_type;pos.Id(ga.pos.Id());
    }

    return *this;
}

const char *Atom::strName() const
{
    switch (my_type) {
        case hydrogen: return "H";
        case carbon: return "C";
        case nitrogen: return "N";
        case oxygen: return "O";
        case sulfur: return "S";
        default: return "Unknown";
    };
}

void Atom::Write()
{
    prf::cout<<"Atom "<<UniqueId()<<"  "<<strName()<<"\t"<<Position()<<'\n';
}

Output & operator<<(Output &op, const Atom &a)
{
    op<<"Atom "<<a.UniqueId()<<"  "<<a.strName()<<"\t"<<a.Position();
    return op;
}
