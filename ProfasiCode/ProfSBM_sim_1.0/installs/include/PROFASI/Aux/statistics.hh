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

#ifndef Statistics_H
#define Statistics_H
#include <cmath>
#include <numeric>
#include <climits>
#include <vector>

namespace prf_utils
{
//! Evaluates the non-normalized correlation &lt;uv&gt;-&lt;u&gt;&lt;v&gt;
//! between to vectors of doubles, u and v.

    class CorrelationEvaluator
    {
    public:
        CorrelationEvaluator() {}

        ~CorrelationEvaluator() {}

        //! Find correlation between vectors v1 and v2
        /**
         * Usage: (for vectors v1 and v2 of double precision values
         * CorrelationEvaluator corr;
         * double correlation= corr(v1,v2);
         */
        double operator()(std::vector<double> & v1,std::vector<double> & v2) {
            double v1sum,v2sum,v1v2sum;
            v1sum=accumulate(v1.begin(),v1.end(),0.0)/v1.size();
            v2sum=accumulate(v2.begin(),v2.end(),0.0)/v1.size();
            v1v2sum=inner_product(v1.begin(),v1.end(),v2.begin(),0.0)/v1.size();
            return (v1v2sum-v1sum*v2sum);
        }
    };

//! Simple utility to find the mean of a vector

    class MeanFinder
    {
    public:
        MeanFinder() {}

        ~MeanFinder() {}

        //! Find mean of the vector v.
        double operator()(std::vector<double> & v) {
            return accumulate(v.begin(),v.end(),0.0)/v.size();
        }
    };

    /*! The StatBox is a little class which is designed to work a bit like
      the old fashioned pocket calculators in the statistical mode. One can reset
      the StatBox, one can put numbers in it, and find mean, standard deviation etc.
      The numbers themselves will not be stored in StatBox. So, correlations can't
      be calculated. But sometimes we are not interested in any more than means and
      variations. To calculate correlations, use CorrelationBox. */

    class StatBox
    {
    private:
        double vmax, vmin;
        double A,Q;
        int numentries,bcor;
    public:
        //! Default constructor
        StatBox()
        {
            reset();
            bcor=1; // Use Bessel's correction by default
        }

        //! Copy constructor, it is possible to copy StatBox objects
        StatBox(const StatBox & stbx)
        {
            vmax=stbx.vmax;
            vmin=stbx.vmin;
            Q=stbx.Q;
            A=stbx.A;
            bcor=stbx.bcor;
            numentries=stbx.numentries;
        }

        //! assignment
        StatBox & operator=(const StatBox &stbx) {
            if (this!=&stbx) {
                vmax=stbx.vmax;
                vmin=stbx.vmin;
                Q=stbx.Q;
                A=stbx.A;
                bcor=stbx.bcor;
                numentries=stbx.numentries;
            }

            return *this;
        }

        ~StatBox() {}

        void use_Bessel_Correction(bool flag)
        {
            bcor=flag?1:0;
        }

        //! Reset StatBox, to start taking statistics afresh.
        void reset()
        {
            A=Q=0;
            numentries=0;
            vmin=1e200;
            vmax=-vmin;
        }

        //! Put a double value x into the StatBox
        void put(double x)
        {
            numentries++;
            Q+=((numentries-1)*(x-A)*(x-A))/numentries;
            A+=(x-A)/numentries;
            vmin=std::min(vmin,x);
            vmax=std::max(vmax,x);
        }

        //! Mean of all values "put" into the StatBox
        double mean() const { return A; }

        //! Maximum value
        double upper() const { return vmax; }

        //! Minimum value
        double lower() const { return vmin; }

        //! Root Mean Square value for all the entries
        double rms() const {
            return sqrt(variance());
        }

        //! Variance
        double variance() const {
            if (numentries<2) return 0;
            else return (Q/(numentries-bcor));
        }

        //! Standard deviation
        double stddev() const { return sqrt(variance()); }

        //! Number of entries
        int numberofentries() const { return numentries; }

        //! Sum of the entries
        double sum() const {return A*numentries; }

        //! Sum of the squares of the entries.
        double sqrsum() const {return Q+A*A;}
    };

//! Small utility class RangeBox.
    /**
     * Throw double precision values into this box and it keeps track of the range
     * of the objects. StatBox also finds the range. So, this class is useful only
     * when StatBox is an overkill.
     */

    class RangeBox
    {
    private:
        double vmax, vmin;
    public:
        RangeBox() {
            reset();
        }

        ~RangeBox() {}

        //!put a value into RangeBox
        void put(double x) {
            vmin=std::min(vmin,x);
            vmax=std::max(vmax,x);
        }

        //! Upper limit of all values
        double upper() { return vmax; }

        //! Lower limit of values
        double lower() { return vmin; }

        //! Reset to start afresh
        void reset() {
            vmin=1e200;
            vmax=-vmin;
        }
    };

//! CorrelationBox finds the correlation coefficient between two data streams
    /**
     * This class finds the correlation cooefficient of two streams
     * of data, without storing individual data points. Similar in usage to the
     * StatBox
     */

    class CorrelationBox
    {
    private:
        double sumx, sumy, sumxy, sumx2, sumy2;
        int numentries;
    public:
        CorrelationBox():sumx(0),sumy(0),sumxy(0),sumx2(0),
                sumy2(0),numentries(0) {}

        ~CorrelationBox() {}

        //! Reset, forget everything, start new.
        void reset() { sumx2=sumy2=sumx=sumy=sumxy=0;numentries=0;}

        //! Put a data snapshot consisting of values from two streams, x and y
        void put(double x, double y) {
            sumx+=x;sumy+=y;
            sumxy+=x*y;
            sumx2+=x*x;
            sumy2+=y*y;
            numentries++;
        }

        //! Take away a point
        /**
         * Caution: Take away can not check if the point was ever put in to
         * the box
         */
        void takeaway(double x, double y) {
            sumx-=x;sumy-=y;
            sumxy-=x*y;
            numentries--;
        }

        //! Get the correlation coefficient
        //! (&lt;xy&gt;-&lt;x&gt;&lt;y&gt;)/(std.dev.(x) * std. dev.(y))
        double correlation() {
            return ((sumxy/numentries)-(sumx/numentries)*(sumy/numentries))/
                   (sqrt((sumx2/numentries-(sumx/numentries)*(sumx/numentries))*
                         (sumy2/numentries-(sumy/numentries)*(sumy/numentries))));
        }
    };
}

#endif
