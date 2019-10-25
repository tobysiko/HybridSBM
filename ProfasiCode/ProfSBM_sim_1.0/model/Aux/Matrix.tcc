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

//#include "Matrix.hh"
template <class T>
Matrix<T>::Matrix() : d1(0), d2(0) {}

template <class T>
Matrix<T>::Matrix(size_t x, size_t y) {
    allocate(x,y);
}
template <class T>
Matrix<T>::~Matrix() {}
template <class T>
Matrix<T>::Matrix(const Matrix &mm) {
    allocate(mm.d1,mm.d2);
    v=mm.v;
}
template <class T>
void Matrix<T>::allocate(size_t x, size_t y) {
    d1=x;d2=y;
    v.resize(x*y);
    v*=0;
}
template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &mm) {
    if (this!=&mm) {
        allocate(mm.d1,mm.d2);
        v=mm.v;
    }
    return *this;
}

template <class T>
T &Matrix<T>::operator()(size_t x, size_t y) {return row(x)[y];}

template <class T>
T Matrix<T>::operator()(size_t x, size_t y) const {return row(x)[y];}

template <class T>
void Matrix<T>::dot(const std::valarray<T> &vec, std::valarray<T> &res)
{
    for (size_t i=0;i<dim1();++i) {
        for (size_t j=0;j<dim2();++j) {
            res[i]+=v[i*d2+j]*vec[j];
        }
    }
}

template <class T>
Matrix<T> & Matrix<T>::operator*=(T d) {(v)*=d;return *this;}
template <class T>
Matrix<T> & Matrix<T>::operator+=(const Matrix<T> &d) {(v)+=((d.v));return *this;}
template <class T>
void Matrix<T>::Write(FILE *fp) {
    fwrite(&d1,sizeof(size_t),1,fp);
    fwrite(&d2,sizeof(size_t),1,fp);
    for(size_t i=0;i<(d1*d2);++i) fwrite(&((v)[i]),sizeof(T),1,fp);
}
template <class T>
void Matrix<T>::Read(FILE *fp) {
    fread(&d1,sizeof(size_t),1,fp);
    fread(&d2,sizeof(size_t),1,fp);
    for(size_t i=0;i<(d1*d2);++i) fread(&((v)[i]),sizeof(T),1,fp);
}
