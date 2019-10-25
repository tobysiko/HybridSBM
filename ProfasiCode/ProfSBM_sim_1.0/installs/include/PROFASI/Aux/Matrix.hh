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

#ifndef Matrix_HH
#define Matrix_HH
#include <valarray>
#include "Slice_iter.hh"
#include <cstdio>

template <class T>

class Matrix
{
protected:
    std::valarray<T> v;
    size_t d1,d2;
public:
    Matrix();
    Matrix(size_t x, size_t y);
    Matrix(const Matrix&);
    Matrix & operator=(const Matrix&);
    Matrix & operator+=(const Matrix&);
    virtual ~Matrix();
    void allocate(size_t x, size_t y);
    inline std::valarray<T> & array(){return v;}

    inline size_t size() const {return d1*d2;}

    inline size_t dim1() const {return d1;}

    inline size_t dim2() const {return d2;}

    inline Slice_iter<T> operator()(size_t i) {return row(i);}

    inline Cslice_iter<T> operator()(size_t i) const {return row(i);}

    inline Slice_iter<T> operator[](size_t i) {return row(i);}

    inline Cslice_iter<T> operator[](size_t i) const {return row(i);}

    inline T get(size_t i, size_t j) {return v[i*d2+j];}

    inline void set(size_t i, size_t j, T x) {v[i*d2+j]=x;}

    inline Slice_iter<T> row(size_t i) {
        return Slice_iter<T>(&v,std::slice(i*d2,d2,1));
    }

    inline Cslice_iter<T> row(size_t i) const {
        return Cslice_iter<T>(&v,std::slice(i*d2,d2,1));
    }

    Slice_iter<T> column(size_t i) {
        return Slice_iter<T>(&v,std::slice(i,d1,d2));
    }

    Cslice_iter<T> column(size_t i) const {
        return Cslice_iter<T>(&v,std::slice(i,d1,d2));
    }

    void dot(const std::valarray<T> &vec, std::valarray<T> &res);

    T & operator()(size_t x, size_t y);
    T operator()(size_t x, size_t y) const;
    Matrix &operator*=(T);
    void Write(FILE *fl);
    void Read(FILE *fl);
};

#include "Matrix.tcc"
#endif
