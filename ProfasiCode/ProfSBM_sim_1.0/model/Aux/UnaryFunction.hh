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

#ifndef UNARYFUNCTION_HH
#define UNARYFUNCTION_HH

class UnaryFunction
{
public:
    UnaryFunction();
    virtual ~UnaryFunction();
    //! Value of the function
    virtual double operator()(double x);
    //! Get a local minimum near a given point
    /**
      The value of the input parameter x0 is shifted to the location of the
      minimum. The return value is the value of the function at the minimum.
      */
    double nearest_minimum(double &x0, double stp);
    void bracket_minimum(double xstart, double step0,
                         double &x0, double &x1, double &x2);
    //! Returns the first derivative at x
    /**
      The base class implements a numerical derivative, using an algorithm
      given in the Numerical Recipes in Fortran. The C++ adaptation here is
      our own.

      The function is virtual, and any derived classes of UnaryFunction are
      free to provide alternative faster derivatives when possible. If the
      function is a closed form algebraic expression in x, it may be best to
      use the corresponding algebraic expression for its derivative.
      */
    virtual double derivative(double x);
    //! Evaluate derivative at x, using st_scale as the starting dx value
    double derivative(double x, double st_scale);
    //! Evaluate derivative at x and return an error estimate in err
    double derivative(double x, double st_scale, double &err);
    //! Set a scale for small (not infinitisimal) changes
    /**
      We have no idea about the relevant ranges of values over which the
      function will change. While evaluating derivatives, we need to make
      small increments in x. But what is small ? With this function, a
      value should be given, specifying a range for small but significant
      change in the function. This is NOT a substitute for the inifinitisimal
      in the definition of the derivative!
      */
    inline void set_small(double x) { h_small=x; }
    //! Get the change in the function at x over a scale h
    /**
      Returns (f(x+h) - f(x-h)). For most applications it is simply subtracting
      the result of two function evaluations. But there are situations where
      the delta can be obtained at a much smaller cost than a single function
      evaluation.
      */
    virtual double delta(double x, double h);
private:
    inline double f(double x) { return (*this)(x); }
    double h_small;
};

#endif // UNARYFUNCTION_HH
