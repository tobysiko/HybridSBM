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

#ifndef FRAGMENTER_HH
#define FRAGMENTER_HH
#include "../Elements/Population.hh"
#include <deque>

class Fragmenter
{
public:
    Fragmenter();
    ~Fragmenter();
    void make_fragments(std::deque<prf_xml::XML_Node *> &frgs);
    void join_fragments(std::deque<prf_xml::XML_Node *> &frgs);
    inline void set_population(prf::Population *popl) {p=popl;}
    inline void default_fragment_size(size_t sz) {batchsize=sz;}
private:
    prf::Population *p;
    size_t batchsize;
};

#endif // FRAGMENTER_HH
