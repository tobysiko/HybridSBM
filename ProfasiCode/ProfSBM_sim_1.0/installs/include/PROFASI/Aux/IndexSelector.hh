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

#ifndef INDEXSELECTOR_HH
#define INDEXSELECTOR_HH
#include <vector>
#include <cstddef>

class IndexSelector
{
public:
    IndexSelector();
    IndexSelector(const IndexSelector &is);
    IndexSelector &operator=(const IndexSelector &is);
    ~IndexSelector();
    void set_range(size_t i1, size_t i2);
    void set_weights(std::vector<double> &wts);
    void init();
    size_t pick(double y);
    inline size_t select(double y) { return pick(y); }
private:
    std::vector<double> p,sp;
    size_t offset;
};

#endif // INDEXSELECTOR_HH
