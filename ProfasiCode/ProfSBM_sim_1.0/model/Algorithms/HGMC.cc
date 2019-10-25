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

#include "HGMC.hh"
#include "../Aux/Constants.hh"
#include "../Aux/fileutils.hh"

using namespace prf;
using namespace UnivConstants;

HGMC::HGMC()
{
    tmin=0.0;
    tmax=1.0;
    set_n_temps(8);
    itmp=0;
    par.option("max_lambda","lmax",1,"(Maximum lambda)");
    par.option("min_lambda","lmin",1,"(Minimum lambda)");
    par.option("num_lambdas","nlmb",1,
                "(Number of lambda values)");
    par.option("lambda_file","lfile",1,"(File with lambda values)");
    //par.disable("temperature");
    tfileinuse=explicitntmps=false;
    lambda_index = 0;

}
/**
  \page gmc_opts Options for multi-temperature simulations
  \li \b --max_temperature or \b -tmax : Maximum temperature. Examples,
  --max_temperature 400 Kelvin or -tmax 400 Kelvin.
  \li \b --min_temperature or \b -tmin : Minimum temperature. Examples,
  --min_temperature 0.58
  \li \b --num_temperatures or \b -ntmp : Number of temperatures. Example:
  -ntmp 16
  \li \b --temperature_file or \b -tfile : Temperature file. Example:
  -tfile T.dat . The format for specifying temperatures in the temperature
  file is explained in the program reference for parallel tempering run
  (\ref partemp_progref). The same format works for simulated tempering and
  simulated annealing.
  */
HGMC::~HGMC() {}

Update *HGMC::perform_update()
{
    return uph.perform_update(itmp);
}

unsigned HGMC::SwitchLambda()
{
    return 0;
}

void HGMC::set_n_temps(size_t i)
{
    if (i<2) {
        prf::cerr<<"HGMC> Refusing to use less than 2 lambda values.\n";
        i=2;
    }
    make_default_lambdas(i);
    uph.set_n_temps(ntmp); // TODO
    //set_n_lambdas(ntmp);
    explicitntmps=true;
}

void HGMC::make_default_lambdas(size_t i)
{
    ntmp=i;
    lambdas.resize(i);
    dndT.resize(i,0);
    if (tmax<tmin) std::swap(tmax, tmin);
    std::cout<<"HGMC::make_default_lambdas(size_t i)\n";
    
    double range = tmax-tmin;
    std::cout<<"ntmp"<<ntmp<<"; tmax="<<tmax<<"; tmin="<<tmin<<"; range="<<range<<"\n";
    for (unsigned c=0;c<ntmp;++c){
    	lambdas[c] = tmin + c*( range/(ntmp-1) );// 1.0/tmax*pow(tmax/tmin,((double)c)/(ntmp-1));
    	std::cout<<"c:"<<c<<" = " <<lambdas[c]<<"\n";
    }
}

void HGMC::set_lambda_range(double t0, double t1)
{
    tmin=t0;
    tmax=t1;
    make_default_lambdas(ntmp);
}

int HGMC::SetLambdas(std::vector<double> tmptemp)
{
    if (tmptemp.size()==0) {
        prf::cerr<<"Read a total of 0 lambda values\n";
        return 0;
    }

    std::sort(tmptemp.begin(),tmptemp.end());

    ntmp=tmptemp.size();
    lambdas.resize(ntmp,0);
    dndT.resize(ntmp,0);

    for (unsigned i=0;i<ntmp;++i) lambdas[i]=tmptemp[i];
    uph.set_n_temps(ntmp);
    //set_n_lambdas(ntmp);
    explicitntmps=true;
    
    std::cout<<"n lambdas="<<ntmp<<"\n";
    return 1;
}

int HGMC::SetLambdas(std::string tempfile)
{
    if (prf_utils::TestFile_r(tempfile.c_str())==0) {
        prf::cerr<<"Reading lambda values from file"<<tempfile<<" failed.\n";
        return 0;
    }

    std::vector<std::string> lines,parts;
    std::string line;
    std::ifstream fin(tempfile.c_str());

    while (getline(fin,line)) lines.push_back(line);
    fin.close();

    //bool invt=false,kelv=false;
    std::vector<double> tmparry;

    for (size_t i=0;i<lines.size();++i) {
        line=lines[i];
        line=prf_utils::trim_str(line);

        if (line.empty()) continue;

        if (line[0]=='#') {
            parts.clear();
            prf_utils::split(line,parts);

            if (parts[0]=="#lambda") {
            	std::cout<<"Found #lambda in file\n";
                //if ((parts[1]=="inverted") && !kelv) invt=true;
                //else if (parts[1]=="Kelvin") kelv=true;
            }
        } else {
            if (prf_utils::is_number(line)) {
                double tmp=strtod(line.c_str(),NULL);
                tmparry.push_back(tmp);
            }
        }
    }

    ntmp=tmparry.size();

    /*if (!invt) {
        for (size_t i=0;i<tmparry.size();++i) {
            if (kelv) tmparry[i]=pru_in_kelvin/tmparry[i];
            else tmparry[i]=1.0/tmparry[i];
        }
    }*/
    Logger(10)<<"Read in "<<ntmp<<" lambdas from "<<tempfile<<"\n";
    
    if (SetLambdas(tmparry)!=0) tfileinuse=true;

    return 1;
}

void HGMC::writeLambdaStat(std::string flnm)
{
    Output ost(flnm.c_str());
    ost <<"# Number of visits to different lambdas ...\n";

    for (unsigned i=0;i<ntmp;++i) ost<<i<<"   "<<dndT[i]<<"\n";

    ost.close();
}

void HGMC::writeLambdas(std::string flnm)
{
    Output fp(flnm.c_str());
    fp<<"# Index lambda\n";

    for (unsigned i=0;i<ntmp;++i) {
        fp<<i<<"\t"<<lambdas[i]<<"\n";
    }

    fp.close();
}

void HGMC::SetLambda(unsigned i)
{

	lambda_index = lambdas[itmp];
   //SetBeta(lambdas[itmp=i%ntmp]);
   dndT[itmp]++;
}

int HGMC::prepare_continue(std::string flnm)
{
    if (prf_utils::TestFile_r(flnm)==0) {
        prf::cerr<<"Reading file "<<flnm<<" failed. Statistics on the "
                <<"number of visits to different lambdas will contain "
                <<"information only from the present run. \n";
        return 0;
    }

    std::ifstream fin(flnm.c_str());
    std::string line;
    getline(fin,line);
    size_t dummy, nvis;
    while (fin>>dummy) {
        fin>>nvis;
        if (dummy<ntmp) dndT[dummy]=nvis;
    }
    fin.close();
    return 1;
}

int HGMC::parseCommand(InstructionString s)
{
    if (s.head()=="lfile" or s.head()=="lambda_file") {
        SetLambdas(s.tail().str());
    } else if (s.head()=="nlmb" or s.head()=="num_lambdas") {
        if (!tfileinuse) set_n_temps(atoi(s.part(1).c_str()));
        std::cout << "nlmb="<<lambdas.size()<<"\n";
    } else if (s.head()=="lmax" or s.head()=="max_lambda") {
        if (!tfileinuse) {
            tmax=prf_utils::stringtodouble(s.tail().str());
            set_lambda_range(tmin,tmax);
            std::cout << "lmax="<<tmax<<"\n";
        }
    } else if (s.head()=="lmin" or s.head()=="min_lambda") {
        if (!tfileinuse) {
            tmin=prf_utils::stringtodouble(s.tail().str());
            if (tmin<=0) tmin=1.0e-7;
            std::cout << "lmin="<<tmin<<"\n";
            set_lambda_range(tmin,tmax);
        }
    }
    //std::cout << "lmin="<<tmin<<"; lmax="<<tmax<<"; nlmb="<<lambdas.size()<<"\n";
    return MC::parseCommand(s);
}

void HGMC::print_setup()
{
    prf::cout<<"\n1 Monte Carlo Cycle = "<<CycleLength()
            <<" Elementary Monte Carlo Steps or Updates\n\n";

    uph.print_setup();

    prf::cout<<"Number of Lambdas = "<<ntmp<<"\n";
}

std::string HGMC::ConfSignature()
{
    std::ostringstream ost;
    ost<<"simulation_method "<<descr<<"\n";
    ost<<"nlmb "<<ntmp<<"\n";
    ost<<"lambda_array ";
    for (size_t i=0;i<lambdas.size();++i) {
        ost<<" "<<lambda(i);
    }
    ost<<"\n";
    return ost.str();
}
