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

#ifndef Slice_iter_HH
#define Slice_iter_HH
#include <valarray>

template<class T> class Slice_iter
{
    std::valarray<T> *v;
    std::slice s;
    size_t curr;
    T & ref(size_t i) const {return (*v)[s.start()+i*s.stride()];}

public:
    Slice_iter(std::valarray<T> *vv, std::slice ss) : v(vv), s(ss), curr(0) {}

    Slice_iter(const Slice_iter<T> &sit) : v(sit.v), s(sit.s), curr(sit.curr) {}

    Slice_iter<T> & operator=(const Slice_iter<T> &sit) {
        if (this!=&sit) {
            v=sit.v;
            s=sit.s;
            curr=sit.curr;
        }

        return *this;
    }

    ~Slice_iter() {}

    Slice_iter end() {
        Slice_iter t=*this;
        t.curr=s.size();
        return t;
    }

    size_t size() const {return s.size();}

    Slice_iter &operator++() {curr++;return *this;}

    Slice_iter operator++(int) {Slice_iter t=*this;curr++;return t;}

    T &operator[](size_t i) {return ref(curr=i);}

    T &operator()(size_t i) {return ref(curr=i);}

    T &operator*() {return ref(curr);}
};

template<class T> bool operator==(const Slice_iter<T> &p, const Slice_iter<T> &q)
{
    return p.curr==q.curr && q.s.stride()==q.s.stride() && p.s.start()==q.s.start();
}

template<class T> bool operator!=(const Slice_iter<T> &p, const Slice_iter<T> &q)
{
    return !(p==q);
}

template<class T> bool operator<(const Slice_iter<T> &p, const Slice_iter<T> &q)
{
    return p.curr<q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start();
}

template<class T> class Cslice_iter
{
    const std::valarray<T> *v;
    std::slice s;
    size_t curr;
    const T & ref(size_t i) const {return (*v)[s.start()+i*s.stride()];}

public:
    Cslice_iter(const std::valarray<T> *vv, std::slice ss) : v(vv), s(ss), curr(0) {}

    ~Cslice_iter() {}

    Cslice_iter end() {
        Cslice_iter t=*this;
        t.curr=s.size();
        return t;
    }

    size_t size() const {return s.size();}

    Cslice_iter &operator++() {curr++;return *this;}

    Cslice_iter operator++(int) {Cslice_iter t=*this;curr++;return t;}

    const T &operator[](size_t i) {return ref(curr=i);}

    const T &operator()(size_t i) {return ref(curr=i);}

    const T &operator*() {return ref(curr);}
};

template<class T> bool operator==(const Cslice_iter<T> &p, const Cslice_iter<T> &q)
{
    return p.curr==q.curr && q.s.stride()==q.s.stride() && p.s.start()==q.s.start();
}

template<class T> bool operator!=(const Cslice_iter<T> &p, const Cslice_iter<T> &q)
{
    return !(p==q);
}

template<class T> bool operator<(const Cslice_iter<T> &p, const Cslice_iter<T> &q)
{
    return p.curr<q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start();
}

#endif
