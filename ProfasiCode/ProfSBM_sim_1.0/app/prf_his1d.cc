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

#include <Aux/fileutils.hh>
#include <Aux/ProgUtils.hh>
#include <Aux/AdaptiveHis.hh>
#include <Aux/profasi_io.hh>
#include <iostream>
#include <string>
using std::string;
using namespace prf;
using namespace prf_utils;

int main(int argc, char *argv[])
{
    prf_utils::ProgArgs par;
    par.option("range","r",2);
    par.option("number_of_bins","nb",1);
    par.option("number_of_blocks","nk",1,"(Number of blocks in the histogram)");
    par.option("adjust_interval","ai",1,"(range adjustment interval)");
    par.option("normalization_mode","nm",1,
        "(possible values: 0, 1 or 2, see below for explanation.)");
    par.option("layout","l",1,"(1 for the old PROFASI 1.1 histogram format)");
    par.option("output_file","o",1);
    par.new_switch("no_adjust","na",false,"(no auto range adjustments)");
    par.analyze(argc,argv);
    if (argc==1 || (!par.option_given("r")) ) {
        prf::cout<<"\nThis program reads data from the standard input and makes "
                <<"a PROFASI style histogram out of it. \n\n";
        prf::cout<<"Typical usage: \n\n";
        prf::cout<<"awk '{print $23;}' datafile | "<<argv[0]
                <<" -nb 100 -r 20 60 -o col23his.dat \n\nor\n\n"
                <<"awk '{print $2,$23;}' datafile | "<<argv[0]
                <<" -nk 8 -nb 100 -r 20 60 -o col23his.dat \n\n";
        prf::cout<<"Options could be any of those described below. "
                <<"The range option \"-r\" is mandatory.\n";
        par.write_available();
        prf::cout<<"For more detailed usage explanation please refer to the "
                <<"documentation.\n\n";
        return 1;
    }
    string ofile="histogram.dat";
    int nbins=200,nadjust=10000,normtype=2,nblocks=1,lyout=2;
    double xmin=0,xmax=1;
    bool fixedhis=par.switch_given("na");
    if (par.option_given("nb")) nbins=atoi(par.option("nb").c_str());
    if (par.option_given("nk")) nblocks=atoi(par.option("nk").c_str());
    if (par.option_given("o")) ofile=par.option("o");
    if (par.option_given("ai")) nadjust=atoi(par.option("ai").c_str());
    if (par.option_given("r")) {
        xmin=strtod(par.option_arr("r",0).c_str(),NULL);
        xmax=strtod(par.option_arr("r",1).c_str(),NULL);
    }
    if (par.option_given("nm")) normtype=atoi(par.option("ntype").c_str());
    if (par.option_given("l")) {
        if (par.option("l")=="1" or par.option("l")=="profasi_his_v1") lyout=1;
    }

    double dummy;
    AdaptiveHis his;
    his.Name(ofile);
    his.Range(xmin,xmax);
    his.Nbins(nbins);
    his.NBlocks(nblocks);
    his.init();
    int idata=1,iblk=0;
    std::cin.tie()->flush();
    if (par.option_given("nk")) {
        while (std::cin>>iblk) {
            std::cin>>dummy;
            his.put(dummy,iblk);
            if (!fixedhis and idata%nadjust==0) his.adjust();
            ++idata;
        }
    } else {
        while (std::cin>>dummy) {
            his.put(dummy,0);
            if (!fixedhis and idata%nadjust==0) his.adjust();
            ++idata;
        }
    }
    if (!fixedhis) his.adjust();
    his.Export(ofile.c_str(),normtype,lyout);
    return 0;
}

/**
  \page prf_his1d Creating PROFASI style histograms with prf_his1d
  You have a data file with hundreds of thousands of lines and tens of columns,
  and you want to look at the distribution of the data in the 13th column.
  prf_his1d is a convenient tool to make such a histogram. You need to use
  UNIX commands "awk" or "cut" to extract the relevant column, and pipe it
  to prf_his1d like this:\n\n

  awk '{print $13;}' the_datafile | prf_his1d -r -30 55 -nb 100 -o output.his \n\n

  The available options are :\n\n

  \li <b>"output_file" or "o"</b>
  \li <b> "range" or "r" </b>: takes two numeric arguments e.g. "-r 1 4.5"
  \li <b>"number_of_bins" or "nb"</b>: takes one integer argument
  \li <b>"normalization_mode" or "nm"</b>: choose a normalisation mode for the histogram.
        possible values: 0, 1 or 2
  \li <b>"layout" or "l"</b>: Use this option to store data in the old PROFASI 1.1
        histogram format. For explanation of data layouts and normalisation modes
        please refer to the documentation of prf_utils::His1D
  \li <b>"no_adjust" or "na"</b>: histogram will not adjust range according to the data
  \li <b>"adjust_interval" or "ai"</b> : takes one integer argument, which gives the number
        of data points to be processed before each range adjustment
  \li <b>"number_of_blocks" or "nk"</b>: takes one integer argument, the number
       of blocks in the histogram (See below for explanation)


  The program simply takes the data from the standard input and puts them in
  a histogram with the specified range (-r -30 55) and number of bins(-nb 100).
  If the data falls outside the specified range, by default, the histogram will
  try to adjust its range to fit the data, without losing information. This is
  a property of PROFASI's self adjusting histograms. With the option "--no_adjust
  or -na" this feature is turned off. Even if you want to use the self adjust
  feature, try to provide at least a sensible range in the range option. This
  is used to set the bin size, which is not changed while adjusting the histogram
  range later.\n\n

  Suppose you have stored some sort of integer identifier in the 2 column,
  that labels the line of data in some way. In PROFASI run time history files,
  you have the temperature index written in the second column. You might be
  interested in the histogram of energy for temperature index 4. You can do
  that exactly as before by using a conditional in "awk".\n\n

  awk '{if ($2==4) print $3;}' n* /rt | prf_his1d -r 10 60 -nb 100 -o etot.his \n\n

  But then you might also need to get the histogram for temperature index 0,1,2... !
  A new file and a new command for every index is unnecessary. PROFASI histograms
  can store, in one file, information about the same kind of data categorized by some
  integer index. The right way to make a single output file containing histograms
  of column 3 partitioned by column 2 is : \n\n

  awk '{print $2,$3;}' n* /rt | prf_his1d -nk 16 -o etot.his -r 10 60 -nb 100 \n\n

  The option "-nk" used in the above example changes how the incoming data stream
  is interpreted. First, the data is now assumed to be in two columns. The first
  column is supposed to be an integer index, whereas the second column is the
  data to be put in the histograms. The value "16" passed to the "-nk" option
  instructs the program to create 16 blocks for the histogram, i.e., that the
  integer index will vary from 0 to 15. The data lines in which the integer index
  is 3 will be stored in the 3rd block of the histogram, index 4 ==> block 4 etc.
  The resulting histogram file will contain many columns. The first column is the
  x bins. The subsequent columns will be the histogram values (frequencies,
  probabilities ...) for the different integer indexes. You can change the
  layout of the histogram to PROFASI 1.1 histogram layout by using the "--layout or -l 1"
  option.

  Although we used the total energy in the above examples, if you are interested
  in that histogram, it is much better to use the histograms generated during the
  runs, like n0/his_Etot etc. This program is intended to be used in situations where
  such a histogram file is not available.


  */
