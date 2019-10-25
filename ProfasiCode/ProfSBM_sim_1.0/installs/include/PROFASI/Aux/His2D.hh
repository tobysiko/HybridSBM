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

#ifndef His2D_HH
#define His2D_HH
#include <cstdio>
#include <string>
#include "profasi_io.hh"
#include <fstream>
#include "Matrix.hh"
#include <vector>

namespace prf_utils
{
    //! A histogram with 2 independent variables
    /**
    * His2D represents histograms of the type dP/dxdy to chart out
    * the probability distribution of a variable z as a function of
    * two variables x and y.
    *
    * Data can be stored in two supported layouts in text files. To
    * describe we use the notation zij=P(xi,yj).
    *
    * The layout means, when the histogram data is written to a file,
    * it is formatted like this:\n \n
    * x0 y0 z00 \n
    * x0 y1 z01 \n
    * x0 y2 z02 \n
    * ... \n
    * x0 yn z0n \n
    * \n
    * x1 y0 z10 \n
    * x1 y1 z11 \n
    * ... \n
    *
    * In the second supported data layout, the data is saved in 3 files. One
    * representing just the x values, { x0, x1,... }, one with only the y
    * values {y1,y2 ...}, and one with the z values :\n\n
    * z00 z01 z02 ... z0n\n
    * z10 z11 z12 ... z1n\n
    * z20 ...\n
    * ...\n
    * where n is (n_ybins -1)
    * This format takes less space on disk. It became the default layout in PROFASI v 1.4.
    *
    *\ingroup utilities
    */

    class His2D
    {
    public:
        His2D();
        His2D(const His2D &);
        ~His2D();
        //! One can assign a histogram to another
        His2D & operator=(const His2D &);
        //! Necessary to call this after specifying number of bins etc.
        void init();
        //! Does not deallocate memory, but clears the collected statistics.
        void reset();
        //! Set range in x and y directions
        inline void Range(double x0,double x1,double x2,double x3)
        {xmin=x0;xmax=x1;ymin=x2;ymax=x3;}

        //! Set only the x range
        inline void XRange(double x0, double x1) {xmin=x0;xmax=x1;}

        //! Set only the y range
        inline void YRange(double x0, double x1) {ymin=x0;ymax=x1;}

        //! Set the number of bins in the X-axis
        inline void NXbins(int v) {nxbin=v;}

        //! Set the number of bins in the Y-axis
        inline void NYbins(int v) {nybin=v;}

        //! Return the current number of bins in the X-axis
        inline int NXbins() const{return nxbin;}

        //! Return the current number of bins in the Y-axis
        inline int NYbins() const{return nybin;}

        //! Return the lower limit of the current X-range
        inline double Xmin() const {return xmin;}

        //! Return the upper limit of the current X-range
        inline double Xmax() const {return xmax;}

        //! Return the lower limit of the current Y-range
        inline double Ymin() const {return ymin;}

        //! Return the upper limit of the current Y-range
        inline double Ymax() const {return ymax;}

        //! X value for the i'th bin on the X-axis
        inline double xval(int i) {return X[i];}

        //! Y value for the i'th bin on the Y-axis
        inline double yval(int i) {return Y[i];}

        //! Z value for the i'th bin along X and j'th bin along Y-axis
        inline double val(int i, int j) {return xy[i][j];}

        //! Put a new datapoint (x,y) into the histogram
        /**
        * We are concerned with keeping track of how often an certain combination
        * (x,y) occurs. This function should be called for every sampling event.
        * Then after a large number of sampling events, the histogram will
        * contain the probabilities (frequencies, if not normalized) for the
        * different (x,y) pairs.
        */
        int put(double x, double y);
        //! Convert the collected frequency data to probabilities
        /**
        * After a call to normalize, the sum of all elements in the histogram
        * is 1, provided all attempted put(x,y) operations succeeded, i.e.,
        * there were no out of range values.
        */
        double normalize();
        //! Save histogram information in a file for plotting
        void Export(const char *filename,int fmt=3);
        //! Read back histogram information from exported text files
        /**
          Import is only implemented for files exported with the default format,
          fmt=3. In this format, there are a few comment lines beginning with
          "#" giving ranges and number of bins. Then there is a matrix M of data
          of (nxbins+1) rows and (nybins+1) columns. The top-left element, with
          index (0,0) is always 0. The elements M[0,1+j], for j in the range
          (0,nybins) are the centres of the ybins for the data. Similarly,
          M[1+i,0] for i in range (0,nxbins) are the centres of xbins. The
          elements M[1+i,1+j] with i in range (0,nxbins) and j in range (0,nybins)
          then represent the probabilities P(X[i],Y[j]).
          */
        int Import(const char *filename);

        //! Save current state of sampling in a binary file to resume counting later
        void save_state(const char *filename);
        //! Retrieve stored sampling state to resume counting
        void read_state(const char *filename);
        //! Add the contents of another histogram.
        His2D & operator+=(const His2D &);
        //! For each x, return Sum z(x,y) | y=0..ymax
        void Projection_X(std::vector<double> &hsx);
        //! For each y, return Sum z(x,y) | x=0..xmax
        void Projection_Y(std::vector<double> &hsy);
    private:
        Matrix<double> xy;
        std::valarray<double> X,Y;
        double xmin, xmax, ymin, ymax,xbin,ybin;
        long nentries,nentries_in;
        int nxbin,nybin;
    };
}

#endif
