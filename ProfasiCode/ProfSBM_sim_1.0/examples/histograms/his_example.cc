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

#include <PROFASI/Aux/fileutils.hh>
#include <PROFASI/Aux/ProgUtils.hh>
#include <PROFASI/Aux/AdaptiveHis.hh>
#include <PROFASI/Aux/profasi_io.hh>
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
    par.option("adjust frequency","na",1,"(range adjustment frequency)");
    par.option("normalization_type","ntype",1,
        "(possible values: 0, 1 or 2, see below for explanation.)");
    par.option("output_file","o",1);
    par.analyze(argc,argv);
    if (argc==1 || (!par.option_given("r")) ) {
        prf::cout<<"This program reads data from the standard input and makes "
                <<"a histogram out of it. \n";
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] input_trajectory_file\n\n";
        prf::cout<<"Options could be any of those described below. "
                <<"The range option \"-r\" is mandatory.\n";
        par.write_available();
        prf::cout<<"\n\nThe range adjustment frequency option (\"-na\") sets "
                <<"the interval at which the histogram attempts to adjust its "
                <<"range according to the data already entered. The default "
                <<"is 10, useful to see the effect if you are typing in the "
                <<"values by hand. If you are piping in the numbers from another "
                <<"program, this option should be given with a value like 1000 "
                <<"to avoid too frequent range adjustments. But don't make it "
                <<"a billion!! In the worst case scenario, this is also the "
                <<"size of out of range list, which is saved until the range "
                <<"adjusts.\n";
        prf::cout<<"\n\nThe normalization type option (\"-ntype\") controls "
            <<"what is printed in the histogram output file. A value of 0 chooses"
            <<" raw frequencies, without normalization. \n"
            <<"The other options normalize to 1 or to the normalization fraction "
            <<"Nf=(1-(n_out_of_range/n_data_entry_attempts)).\n"
            <<"A value of 1 for ntype normalizes the bins, such that the  "
            <<"values for the different bins add up to the normalization "
            <<"fraction. A value of 2 means, sum of "
            <<"bin values times bin widths is the normalization fraction.\n";
        return 1;
    }
    string ofile="histogram.dat";
    int nbins=200,nadjust=10,normtype=2;
    double xmin=0,xmax=1;
    if (par.option_given("nb")) nbins=atoi(par.option("nb").c_str());
    if (par.option_given("o")) ofile=par.option("o");
    if (par.option_given("na")) nadjust=atoi(par.option("na").c_str());
    if (par.option_given("r")) {
        xmin=strtod(par.option_arr("r",0).c_str(),NULL);
        xmax=strtod(par.option_arr("r",1).c_str(),NULL);
    }
    if (par.option_given("ntype")) normtype=atoi(par.option("ntype").c_str());
    double dummy;

    AdaptiveHis his;
    his.Range(xmin,xmax);
    his.Nbins(nbins);
    his.init();
    int idata=1;
    std::cin.tie()->flush();
    while (std::cin>>dummy) {
        his.put(dummy);
        if (idata++%nadjust==0) his.adjust();
    }
    his.adjust();
    his.Export(ofile.c_str(),normtype);
    return 0;
}

