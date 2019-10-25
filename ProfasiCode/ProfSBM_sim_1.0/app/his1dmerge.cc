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

#include <Aux/AdaptiveHis.hh>
#include <vector>
#include <algorithm>
#include <Aux/fileutils.hh>
#include <Aux/ProgUtils.hh>
#include <Aux/profasi_io.hh>

using std::vector;
using std::find;
using std::string;

using namespace prf;

using namespace prf_utils;

int main(int argc, char *argv[])
{
    prf_utils::ProgArgs par;
    par.option("output_file_name","o",1);
    par.option("normalization_mode","nm",1,
               "(normalization mode of the output histogram)");
    par.option("layout","l",1,"(data layout in the output file)");
    par.option("max_load","mx",1,"(Maximum number of histograms to load simultaneously)");
    par.new_switch("verbose","v",false,"(More verbose logs)");
    par.disable("st");
    par.analyze(argc,argv);

    if (argc==1) {
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [OPTIONS] input_file(s)\n\n";
        prf::cout<<"Options could be any of ...\n";
        par.write_available();
        prf::cout<<"The program prints the names of the files it is working "
                <<"on, if verbose mode is set.\n";
        prf::cout<<"For explanation on the normalization modes and data layout "
                <<"in ProFASi histograms, please see the documentation of class "
                <<"His1D. The valid valued for normalization mode for option "
                <<"\"-nm\" are 0, 1 or 2. For option \"-l\" the valid values are "
                <<"1 or 2. 1 means format profasi_his_v1, 2 means format "
                <<"profasi_his_v2. Unrecognized values lead to default values "
                <<"2 for normalization mode and profasi_his_v2 for layout.\n";
        return 1;
    }

    bool verbose_mode=false;
    int nrmmd=2,lyout=2;

    string ofile="result";
    vector<string> filename;

    if (par.switch_given("v")) verbose_mode=true;
    if (par.option_given("o")) ofile=par.option("o");
    if (par.option_given("nm")) nrmmd=atoi(par.option("nm").c_str());
    if (nrmmd!=1 && nrmmd!=0) nrmmd=2;
    if (par.option_given("l")) {
        if (par.option("l")=="1" or par.option("l")=="profasi_his_v1") lyout=1;
    }


    for (int i=0;i< par.n_spare_args();++i) {
        string fcandidate=par.spare_args(i);
        if (TestFile(fcandidate)) filename.push_back(fcandidate);
    }

    if (filename.empty()) {
        prf::cerr<<"no input files\n";
        return 1;
    }

    unsigned int mxloads=filename.size();
    if (par.option_given("mx")) mxloads=strtoul(par.option("mx").c_str(),NULL,10);
    if (mxloads>filename.size()) mxloads=filename.size();

    vector<AdaptiveHis> his(mxloads);
    AdaptiveHis sumhis;
    sumhis.Name(ofile);

    size_t i=0;
    while (i<filename.size()) {
        if (verbose_mode) prf::cout<<"file : "<<filename[i]<<"\n";

        his[i%mxloads].Import(filename[i].c_str());
        if ((i+1)%mxloads==0 or (i+1)==filename.size()) {
            AdaptiveHis partialsum=add(his);
            if ((i+1)==mxloads) sumhis=partialsum;
            else sumhis+=partialsum;
        }
//        his[i].unnormalize();
        ++i;
    }

//    sumhis.normalize();
    sumhis.Export(ofile.c_str(),nrmmd,lyout);
    return 0;
}

/**
\page his1dmerge Merging and converting ProFASi histograms with his1dmerge
The program "his1dmerge" reads one or more histogram files, combines
all read histograms and creates a single histogram file out of their data.
His1D or AdaptiveHis objects store normalization and range information.
These can be used to construct properly normalized total histograms.
The data layout and normalization mode of the output histograms can be
chosen here. This means, this program also functions as a converter
for different histogram forms. The layout of the input histogram will
be inferred from what is in there. Only the layout of the output histogram
can be passed as an argument. Histograms of different layout can be
merged without problems.<br>
\sa prf_utils::His1D, prf_utils::AdaptiveHis

Syntax:<br><br>

his1dmerge.ex [OPTIONS] input_file(s)<br><br>

Options could be any of ...<br>
<ul>
<li><b>-o  or --output_file_name </b> 1 argument, the filename</li>
<li><b>-nm  or --normalization_mode </b> 1 argument, the desired normalization
mode in the output file</li>
<li><b>-l  or --layout </b> 1 argument, the desired data layout in
the output file</li>
<li><b>-v  or --verbose </b>without any arguments. The program prints the
names of the files it is working on, if verbose mode is set. </li>
</ul>
Example:<br><br>

his1dmerge -v -o test.his n* /his_Etot <br><br>

The above will create a combined histogram from all files that match the
pattern n* /his_Etot. Note that His1D objects are intrinsically <i>blocked</i>,
meaning each histogram can contain multiple blocks of data, corresponding to,
for instance, different temperatures.<br><br>
his1dmerge oldhis.dat -nm 2 -l 2 -o newhis.dat <br><br>

The above will convert an old ProFASi histogram to a new one with normalization
mode 2 and layout 2.

For explanation of available normalization modes and different data layouts
please see the prf_utils::His1D class documentation.

*/
