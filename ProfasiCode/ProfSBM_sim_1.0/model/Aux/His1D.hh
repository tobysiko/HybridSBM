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

#ifndef His1D_HH
#define His1D_HH
#include "Matrix.hh"
#include <vector>
#include <string>

namespace prf_utils
{
    //! The default histogram utility in ProFASi is based on His1D
    /**
     * This class can be used to construct histograms of one variable in ProFASi.
     * In practice, the derived class AdaptiveHis, which can adjust its own
     * range is used. But most of the functionality of the histograms is
     * implemented in His1D, and is described here.
     * Typical use:
     \verbatim
       His1D his;
       his.Range(0,25);
       his.NBins(100);
       his.init();
       while (you are doing something useful) {
          do_some_great_calculations();
          x=a_very_important_number();
          his.put(x);
       }
       his.Export("his_avim.dat");
     \endverbatim
     * When you "put" lots of numbers into a His1D object as above, it measures
     * how often the values you put were in different bins. The bins are
     * calculated between the range you provide and the number of bins you set.
     * Upon "Export()", you get a file with what the histogram has measured.
     * Just how the measured histograms are reported in the output file varies,
     * and needs a bit more explanation.
     *
     * <h3> Normalization </h3>
     * A histogram may or may not be normalized. If it is not normalized, the y
     * values for each x bin is simply the frequency of data points in the range
     * (x-0.5*dx, x+0.5*dx). There is a function normalize() that can ask the
     * histogram to normalize its data, so that the sum of the y values for all
     * x bins is 1. The normalization can be undone with the unnormalize()
     * function. But an explicit call to normalize() is normally not necessary.
     * When the histogram is saved, the data written to the file is normalized,
     * without changing the internal frequency counts. This way, a histogram can
     * save normalized data at intermediate stages during a run, without messing
     * up the frequencies. So, you don't normally need to call normalize().
     *
     * Sometimes, the absolute observed frequencies rather than normalized data
     * is of interest. So, it should be possible to save them. Therefore, the
     * function "Export()" can be given a normalization mode parameter.
     *
     * \c his.Export("his_avim.dat",0);
     *
     * The raw frequencies are then saved with no normalization. Two additional
     * modes of normalization are supported. Normalization mode 1 means that
     * the sum of the values at all bins is 1. In normalization mode 2, sum of
     * the values at each bin multiplied by the bin sizes is 1. In mode 1,
     * two histograms of the same data, with different bin sizes, will
     * superficially look different when plotted. The one with bigger bins will
     * end up with larger amplitudes. In mode 2, the height of the histogram
     * plots will be unchanged when bin size is changed.
     *
     * <h3>Histogram data blocks</h3>
     * Another important aspect of His1D objects is that they are intrinsically
     * indexed. Meaning, if you need 8 histograms for the same quantity
     * corresponding to 8 values of some state designator, you only need
     * one His1D. But you give it 8 blocks with the NBlocks(int) function.
     * As an example, we take a parallel tempering run with 64 temperatures.
     * Histograms of energy at every temperature is interesting. But then we
     * end up with 64 histograms for total energy, and 64 again for every
     * other quantity we wish to measure. In situations like this, you create
     * a His1D object with many blocks, 64 in this case. Data for
     * different blocks are written into the same output file, and of course,
     * each block is separately normalized. In principle, this is just a way
     * to combine related histograms into one file for the output, and one
     * object for less cumbersome handling in the program.
     *
     * <h3>Data layout in ProFASi histogram files</h3>
     *
     * There are two supported layouts for the histogram data. One can choose
     * the desired layout by passing an additional parameter to "Export()".
     *
     * \c his.Export("his_avim.dat", 2,2);
     *
     * The layout referred to as "profasi_his_v1" is the layout familiar to
     * ProFASi users from older versions of the program. The data blocks
     * are organized as:
     \verbatim
       b0 x0 y00
       b0 x1 y01
       b0 x2 y02
       b0 x3 ...
       ...
       blank line
       b1 x0 y10
       b1 x1 y11
       etc.
     \endverbatim
     * In the above, b0, b1 etc were integer block identifiers. The different
     * blocks could represent, for instance, different temperatures. Plotting
     * this kind of a histogram file in GNUPLOT is easy:
     *
     * \c gnuplot> p "myhis.dat" u 2:3 w l
     *
     * This gives one line for each block. But this format is not very economic.
     * The X values are repeated in every block. The numbers b0, b1 etc are
     * written many times. The second supported format addresses these issues.
     * The layout referred to as "profasi_his_v2" is like this:
     \verbatim
       x0 y00 y10 y20 y30 ...
       x1 y01 y11 y21 y31 ...
       x2 y02 y12 y22 y31 ...
       etc.
     \endverbatim
     * The first column has the X values. The subsequent columns have the Y
     * data for blocks 0, 1, 2 ... The block numbers are never written. The
     * X values are written once. To plot it in GNUPLOT write:
     *
     * \c gnuplot> p "myhis.dat" u 1:someblock w l
     *
     * This will plot only one line at a time, depending on what block you
     * choose.
     *
     * If the histogram has only a single block, the block id is always 0
     * and is not written even in the layout profasi_his_v1. In this
     * situation, the two layouts are equivalent. Histograms generated
     * with ProFASi versions up to 1.2 were always in layout profasi_his_v1.
     * In version 1.5, the default layout is changed to profasi_his_v2.
     * The program his1dmerge (See \ref his1dmerge) can be used to change
     * histogram data layout as well as normalization modes.
     *
     * \sa his1dmerge
     * \ingroup utilities
     */

    class His1D
    {
    public:
        //! Default constructor
        His1D();
        //! Create His1D with nbl blocks (create nbl histograms)
        His1D(int nbl);
        //! Construct with xmin, xmax, number of bins and number of blocks
        His1D(double xmn,double xmx,int npnts,int numblocks=1);
        //! Construct with xmin, xmax, bin size and number of blocks
        /**
          The given range (xmn, xmx) is symmetrically extended to accommodate
          an integral number of bins of the requested size, if needed.
          */
        His1D(double xmn, double xmx, double bnsz, int numblocks=1);
        //! Copy constructor (copies data, so do not initialize after this!)
        His1D(const His1D &);
        virtual ~His1D();
        //! Assignment operator (copies data, do not initialize after this!)
        His1D &operator=(const His1D &);
        //! Give it a name
        inline void Name(std::string nm) {nam=nm;}
        inline std::string Name() const {return nam;}
        //! Initializes using info about range etc
        void init();
        //! Reset all data (init calls this)
        void reset();
        //! Make it a histogram of n blocks
        inline void NBlocks(int n) {nblk=n;}

        //! Return the number of blocks
        inline int NBlocks() const {return nblk;}

        //! number of entries in block i
        inline long n_entries(int i) const {return nentries[i];}

        //! number of entries in block i in range
        inline long n_entries_in_range(int i) const {return nentries_in[i];}

        //! Set range
        inline void Range(double x0,double x1) {xmin=x0;xmax=x1;}

        //! Set number of bins
        inline void Nbins(int v) {nbin=v;}

        //! Get number of bins
        inline int Nbins() const{return nbin;}

        //! Set bin size
        /**
          Create bins of a given size. The existing histogram range is extended
          to accommmodate an integral number of bins.
          */
        void set_bin_size(double sz);

        //! Get xmin
        inline double Xmin() const{return xmin;}

        //! Get xmax
        inline double Xmax() const{return xmax;}

        //! Get bin size
        inline double Xbin() const {return xbin;}

        //! x value for the middle of i'th bin
        inline double xval(int i) {return xmin+(i+0.5)*xbin;}

        //! y value for the middle of i'th bin
        inline double yval(int iblk,int i) {return xy[iblk][i];}

        //! put value x into the iblk block
        int put(double x, int iblk=0);
        //! put n indentical values at once
        int nput(double howmanytimes, double x, int iblk);
        //! normalize histogram so that each block sums to 1
        double normalize();
        //! unnormalize histogram so that each block sums to its occupancy
        double unnormalize();
        //! Export data in text format to a given file
        /**
        * Choose between different normalization modes and data layout options
        * as described above.
        */
        virtual void Export(const char *filename, int normmode=2, int datlayout=2);
        inline void Export() {Export(nam.c_str(),2,2);}
        //! Import histogram data written in the format of Export function
        int Import(const char *filename);
        //! Add information from another given histogram
        virtual His1D & operator+=(His1D &);
        void import_data(Matrix<double> &mtx);
        inline Matrix<double> &data() {return xy;}

    protected:
        Matrix<double> xy;
        double xmin, xmax, xbin;
        std::vector<double> nentries,nentries_in;
        int nbin,nblk;
        std::string nam;
    };
}

#endif
