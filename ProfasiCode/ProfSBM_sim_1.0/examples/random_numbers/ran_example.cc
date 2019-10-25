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
#include <PROFASI/Aux/RandomNumber.hh>
#include <PROFASI/Aux/Constants.hh>
#include <iostream>
#include <string>
using std::string;
using namespace prf;

class Gaussian : public RandomNumber {
private:
    double x,y,m,d;
    unsigned long req;
public:
    Gaussian();
    ~Gaussian();
    double shoot();
    inline void set_pars(double mn, double st) {m=mn;d=st;}
};
Gaussian::Gaussian() : RandomNumber() {
    RandomizeState();
    req=0;m=x=y=0; d=1; 
}
Gaussian::~Gaussian() {}

double Gaussian::shoot() 
{
    double r,phi;
    if (req%2==0) {
        r=sqrt(2)*d*sqrt(-log(RandomNumber::shoot()));
        phi=UnivConstants::twoPi*RandomNumber::shoot();
        x=m+r*cos(phi);
        y=m+r*sin(phi);
        return x;
    } else return y;
}


int main(int argc, char *argv[]) 
{
    prf_utils::ProgArgs par;
    par.option("mean","m",1);
    par.option("standard_deviation","d",1);
    par.option("how_many","n",1);
    par.analyze(argc,argv);
    if (!(par.option_given("m") || par.option_given("d"))) {
        std::cout<<"This program prints a series of random numbers distributed "
                <<"as a Gaussian of given mean and standard deviation. \n";
	std::cout<<"Usage: \n\n"; 
	std::cout<<argv[0]<<" [OPTIONS] \n";
	std::cout<<"Options could be...\n";
	par.write_available();
	return 1;
    }
    int n=10;
    prf::Logger::verbosity=1;
    double mean=0,stdev=1;
    if (par.option_given("n")) n=atoi(par.option("n").c_str());
    if (par.option_given("m")) mean=strtod(par.option("m").c_str(),NULL);
    if (par.option_given("d")) stdev=strtod(par.option("d").c_str(),NULL);
    Gaussian ran;
    ran.set_pars(mean,stdev);
    for (int i=0;i<n;++i) std::cout<<(ran.shoot())<<"\n";
    return 0;
}

