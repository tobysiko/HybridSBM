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

#include <PROFASI/Aux/ProgUtils.hh>
#include <PROFASI/Aux/fileutils.hh>
#include <PROFASI/Elements/Population.hh>
using namespace std;

using namespace prf;

using namespace prf_utils;

int main ( int argc, char *argv[] )
{
    ProgArgs par;
    par.option ( "log_level","ll",1 );
    par.option ( "output_file","o",1);
    par.analyze ( argc,argv );

    if ( argc==1 ) {
        std::cout<<"Usage: \n\n";
        std::cout<<argv[0]<<" [parameters] sequence\n";
        std::cout<<"parameters could be any number and combination of ...\n";
        par.write_available();
        std::cout<<"Examples: \n\n";
        std::cout<<argv[0]<<" \"ALA TRP ALA\"\n";
        return 0;
    }
    Groups::initGroups();
    AminoAcid::initCommon();
    int loglev=10;
    if (par.option_given("ll")) loglev=atoi(par.option("ll").c_str());
    prf::Logger::verbosity=loglev;
    double boxl=100.0;
    AtomCoordinates::SetBox(boxl);
    string seq="ALA",ofile="new.pdb";

    if ( par.n_spare_args() >=1 ) seq=par.spare_args ( 0 ) ;
    if (par.option_given("o")) ofile=par.option("o");

    Population p;
    p.AddProtein(seq,1);
    p.Init();
    p.InitCoord("stretched");

    FILE *fp=fopen(ofile.c_str(),"w");
    p.WritePDB(fp);
    fclose(fp);
    return 0;
}

