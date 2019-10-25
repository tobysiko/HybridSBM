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

#include <Aux/His2D.hh>
#include <vector>
#include <algorithm>
#include <Aux/fileutils.hh>
#include <Aux/ProgUtils.hh>

using std::vector;
using std::find;
using std::string;

using namespace prf_utils;

int main(int argc, char *argv[])
{
    ProgArgs par;
    par.option("output_file_name","o",1);
    par.new_switch("verbose", "v",false);
    par.option("layout","l",1,"(choose data layout in the output file)");
    par.new_switch("read_binary","b",false,"(Whether you want to read in binary data)");
    par.analyze(argc,argv);

    if (argc==1) {
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] input_file_name(s)\n\n";
        prf::cout<<"Options could be any of the following ...\n";
        par.write_available();
        prf::cout<<"\nThe -l option for the data layout could be either 1 or 2.\n"
                <<"Please refer to the documentation for explanation of the data layouts.\n\n";
        prf::cout<<"Example:\n\n"
                <<argv[0]<<" n*/BackboneRMSD_Etot_4 -o test.his --verbose\n\n";
        return 1;
    }

    vector<string> filename;
    bool verbose_mode=false,bininput=false;
    string ofile="result";
    int fmt=3;

    if (par.switch_given("v")) verbose_mode=true;
    if (par.switch_given("b")) bininput=true;

    if (par.option_given("o")) ofile=par.option("o");

    if (par.option_given("l")) fmt=atoi(par.option("l").c_str())&1;

    for (int i=0;i<par.n_spare_args();++i) {
        string fcandidate=par.spare_args(i);

        if (TestFile(fcandidate)) filename.push_back(fcandidate);
    }

    if (filename.empty()) {
        prf::cerr<<"no input files\n";
        return 1;
    }

    His2D contour;

    if (verbose_mode) prf::cout<<"file : "<<filename[0]<<"\n";

    if (bininput) contour.read_state(filename[0].c_str());
    else contour.Import(filename[0].c_str());

    for (size_t i=1;i<filename.size();++i) {
        if (verbose_mode) prf::cout<<"file : "<<filename[i]<<"\n";

        His2D nxt;
        if (bininput) nxt.read_state(filename[i].c_str());
        else nxt.Import(filename[i].c_str());
        contour+=nxt;
    }

    contour.normalize();
    contour.Export(ofile.c_str(),fmt);
    return 0;
}

/**
\page his2drender Merging binary data from His2D objects
The program his2drender.ex creates a plotable file out of the
binary data saved by a His2D object in PROFASI. Like his1dmerge, it can
merge information from several histograms and normalise the result, if
multiple files are passed. Two data layouts are supported at the moment.
For details on the data layout, please refer to the documentation of
prf_utils::His2D. The default layout is 2. <br>

Syntax:<br><br>

his2drender [OPTIONS] input_file(s)<br><br>
Options could be any of ...<br>
<ul>
    <li><b>-o  or --output_file_name </b> 1 argument, the filename</li>
    <li><b>-l or --layout</b> 1 or 2</li>
    <li><b>-v  or --verbose </b>without any arguments. The program prints the
    names of the files it is working on, if verbose mode is set. </li>

</ul>
Example:<br><br>

his2drender n* /BackboneRMSD_Etot_4 -o test.his --verbose<br>
<br>
The above will create a single 3d plot file out of all files that match the
pattern n* /BackboneRMSD_Etot_4 and save the result to the specified
output.

*/
