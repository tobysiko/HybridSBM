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
#include <Aux/His2D.hh>
#include <Aux/profasi_io.hh>
#include <iostream>
#include <string>
using std::string;
using namespace prf;
using namespace prf_utils;

int main(int argc, char *argv[])
{
    prf_utils::ProgArgs par;
    par.option("range1","r1",2,"(range for the first axis)");
    par.option("range2","r2",2,"(range for the second axis)");
    par.option("number_of_bins1","nb1",1,"(number of bins along the first axis)");
    par.option("number_of_bins2","nb2",1,"(number of bins along the second axis)");
    par.option("output_file","o",1);
    par.option("layout","l",1,"(Data layout in the output file, possible values 1,2  or 3)");
    par.analyze(argc,argv);
    if (argc==1 || (!par.option_given("r1")) || (!par.option_given("r2")) ) {
        prf::cout<<"This program reads data from the standard input and makes "
                <<"a PROFASI style 2d histogram out of it. A 2d histogram is "
                <<"a histogram with two independent variables. Each bin represents "
                <<"a rectangle centred around a point (x,y) with a bin width of "
                <<"(dx,dy) along the x and y axis respectively. The histogram counts "
                <<"how many times a data point fell in each bin of a 2d grid of bins. "
                <<"\n\nFor this program, the data is expected to be in a two column "
                <<"format representing the x and y values. \n\n";
        prf::cout<<"Typical use: \n\n";
        prf::cout<<"awk \'{print $3,$13;}' n*/rt | "
                <<argv[0]<<" -r1 25 200 -nb1 100 -r2 0 25 -nb2 100\n\n";
        prf::cout<<"The above should take columns 3 and 13 from all files n0/rt, n1/rt "
                <<"etc and make a 2d histogram with ranges (25,200) and (0,25) along "
                <<"the x and y axes, and 100 bins along each axis. \n\n";
        prf::cout<<"\n\nRemember that PROFASI 2d histograms do not self adjust the ranges. "
                <<"Any out of range data point is simply ignored.\n";
        prf::cout<<"The \"-l\" option can be used to choose between one of the two supported "
                <<"data layouts for this kind of histogram. \n\n "
                <<"Please refer to the documentation for further details.\n";
        par.write_available();
        prf::cout<<"Options r1 as well as r2 are mandatory.\n";
        return 1;
    }
    string ofile="histogram.dat";
    int nbins1=50,nbins2=50;
    double xmin=0,xmax=1,ymin=0,ymax=1;
    if (par.option_given("nb1")) nbins1=atoi(par.option("nb1").c_str());
    if (par.option_given("nb2")) nbins2=atoi(par.option("nb2").c_str());
    if (par.option_given("o")) ofile=par.option("o");
    if (par.option_given("r1")) {
        xmin=strtod(par.option_arr("r1",0).c_str(),NULL);
        xmax=strtod(par.option_arr("r1",1).c_str(),NULL);
    }
    if (par.option_given("r2")) {
        ymin=strtod(par.option_arr("r2",0).c_str(),NULL);
        ymax=strtod(par.option_arr("r2",1).c_str(),NULL);
    }
    int fmt=3;
    if (par.option_given("l")) fmt=(atoi(par.option("l").c_str())&1);

    His2D his;
    his.XRange(xmin,xmax);
    his.YRange(ymin,ymax);
    his.NXbins(nbins1);
    his.NYbins(nbins2);
    his.init();
    double x,y;
    std::cin.tie()->flush();
    while (std::cin>>x) {
        std::cin>>y;
        his.put(x,y);
    }
    his.normalize();
    his.Export(ofile.c_str(),fmt);
    return 0;
}

/**
  \page prf_his2d Creating PROFASI style 2-d histograms with prf_his2d
  This application creates a density plots P(x,y) out of a stream of
  (x,y) values. The data is assumed to be a white space separated stream
  of numbers with no punctuation. You need to use UNIX commands "awk" or
  "cut" to extract the relevant columns of data if the data file contains
  many columns, and pipe it to prf_his2d like this:\n\n

  awk '{print $15,$3;}' the_datafile | prf_his2d -r1 0 25 -nb1 100 -r2 10 200 -nb2 100 -o output.his \n\n

  The available options are :\n\n

  \li <b>"output_file" or "o"</b>
  \li <b> "range1" or "r1" </b>: range in X axis
  \li <b>"number_of_bins1" or "nb1"</b>: number of bins in X axis
  \li <b> "range2" or "r2" </b>: range in Y axis
  \li <b>"number_of_bins2" or "nb2"</b>: number of bins in Y axis
        possible values: 0, 1 or 2
  \li <b>"layout" or "l"</b>: Use this option to control the layout of data in the output file.
         For explanation of supported data layouts please refer to the documentation of
         prf_utils::His2D

  */

