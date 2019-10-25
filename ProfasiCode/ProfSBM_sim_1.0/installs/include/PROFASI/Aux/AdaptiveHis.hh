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

#ifndef AdaptiveHis_HH
#define AdaptiveHis_HH
#include "His1D.hh"
#include <list>

namespace prf_utils
{
    //! A histogram that can adjust its own range according to the data
    /**
    * This is a minor modification of the His1D class that frees the user
    * from estimating a good range for a histogram before starting to
    * fill it. When one starts a Monte Carlo run with a new protein, or
    * in any other application where a histogram may be required, one has
    * often no idea where the values of one measurement might lie. Before
    * version 1.1 of PROFASI, one had to first make a trial run to get
    * a good feeling for the true range of the data, and then start new
    * production runs where the correct histogram ranges were specified.
    * To a large extent, this will now be unnecessary.
    *
    * This histogram follows the data, wherever it is. You just have to
    * declare an AdaptiveHis, fill it with data, ask it to "adjust()"
    * once in a while, and Export() it to a file. If you then plot the file,
    * it will have a very reasonable range: not too many empty bins on left
    * and right. Not too many missed data points.
    *
    * This does not mean that one does not need to think about the size of
    * the values put into the histograms at all. The adjust() function
    * does not change the size of the bins used for the histogram. That's
    * how it works! If there are many points outside the current range,
    * the class remembers those missed points. When statistics is collected
    * we pretend that there were an infinite number of bins of the size
    * set at initialization. We only choose to do the book keeping on a
    * finite range of those bins. If there are points outside our currently
    * tracked bins, we can add a few bins to the left or right to
    * accommodate them, without affecting the collected data at all. If
    * we were to change the bin size during an adjustment, that would
    * interfere with the data collected before adjust() was called.
    *
    * So, an initial guess for the size of the data is useful. More precisely
    * a good initial estimate of the size of the bins is useful. If the
    * minimum or maximum ranges are wrong, this class will take care of it.
    * If your data values range from 3000 to 10000 and you initialize your
    * AdaptiveHis to a range 0 to 1 with 100 bins, The adjust function
    * will result in a huge histogram. It will have a range 3000 to 10000,
    * with 700000 bins. But if you initialize it to 2000 to 4000 with 50
    * bins, you will be fine. It will once again find the correct range,
    * but will have less than 200 bins.
    *
    * This class was introduced in version 1.1. It is possible that in the
    * future, its functionality will be absorbed in His1D.
    *
    * \ingroup utilities
    */

    class AdaptiveHis : public His1D
    {
    public:
        //! Default constructor
        AdaptiveHis();
        //! Copy constructor
        AdaptiveHis(const AdaptiveHis &);
        ~AdaptiveHis();
        //! Assignment operator
        /**
        * It copies data, do not initialize after this!
        */
        AdaptiveHis &operator=(const AdaptiveHis &);
        //! Initializes using info about range etc
        void init();
        void reset();
        //@{
        /** Disable/enable range tracking features.
        * One can temporarily disable the "adaptive" qualities of the
        * histogram. This is intended for use, if it is known that the
        * incomming data for a certain stage of the program can contain
        * non-sensical values which should have no bearing on the range. When
        * "adjustability" is disabled, the histogram forgets new out of
        * range values, until it is re-enabled.
        */
        inline void disable_adjust() {adj_enabled=false;}

        inline void enable_adjust() {adj_enabled=true;}

        //@}
        //! Adjust range to accommodate data, keeping bin size fixed.
        /**
        * Appropriate range for the data is found by examining out of
        * range points, and the occupancy of currently used bins. We
        * add bins, only if we can fill them. The fundamental reason for
        * the existance of any kind of histograms is that one does not
        * wish to save each and every data point. "Similar" data points
        * are groupped, or binned together. Now if there is one data point
        * outside the current range, such that we would need add 100
        * bins to the right to reach that datapoint, the use of the histogram
        * itself loses its meaning, if we do that. It is more economical
        * to save that data point than to create 100 more bins and remember
        * the frequency of each. Therefore, even after repeated calls to
        * adjust, there might be one remaining out-of-range point that
        * this class simply refuses to cover. When the histogram is saved,
        * such points, if any, are saved in a separate file, and you can
        * deal with them if you like.
        *
        * Return value is non-zero if the range really changes.
        */
        int adjust();
        //! Put a value into the histogram
        /**
        * The value x is put into the histogram if it fits in the range.
        * If not, it is stored in the out-of-range list. The adjust function
        * deals with these out of range points and may put them into bins
        * when the range is appropriately extended.
        */
        int put(double x, int i=0);
        int nput(double howmanytimes, double x, int iblk);
        inline std::list<std::pair<double,int> > &out_of_range_list() {return oor;}

        int check_compatibility(const AdaptiveHis &hs);
        friend AdaptiveHis add(std::vector<AdaptiveHis> &hsv);
        friend AdaptiveHis add(std::vector<AdaptiveHis*> hsv);

        AdaptiveHis &operator+=(AdaptiveHis &hs);

        //! Save histogram and out of range points to files
        /**
        * The histogram data is saved as in class His1D. The parameters
        * normmode and lyout are simply passed down to the base class function.
        * But this class also saves the out-of-range points not covered by
        * the final range(those the function adjust() refuses to include), in
        * a second file with the same name as the histogram file, but with
        * an extension ".out_of_range" at the end.
        */
        void Export(const char *filename, int normmode=2, int lyout=2);
    private:
        bool is_good_bin_boundary(double somept);
        int adapt_bins(int newmin, int newmax);
        int sort_oorpoints();
        int find_new_bin_range(int &newmin, int &newmax);
        std::list<std::pair<double,int> > oor; //oor = out of range
        bool adj_enabled;
    };
}

#endif
