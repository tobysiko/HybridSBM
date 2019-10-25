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

#include "GMC.hh"
#include "../Aux/Constants.hh"
#include "../Aux/fileutils.hh"

using namespace prf;
using namespace UnivConstants;

GMC::GMC()
{
    tmin=200/UnivConstants::pru_in_kelvin;
    tmax=1000/UnivConstants::pru_in_kelvin;
    set_n_temps(8);
    itmp=0;
    par.option("max_temperature","tmax",1,"(Maximum temperature)");
    par.option("min_temperature","tmin",1,"(Minimum temperature)");
    par.option("num_temperatures","ntmp",1,
                "(Number of temperatures)");
    par.option("temperature_file","tfile",1,"(File with temperatures)");
    par.disable("temperature");
    tfileinuse=explicitntmps=false;
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
GMC::~GMC() {}

Update *GMC::perform_update()
{
    return uph.perform_update(itmp);
}

unsigned GMC::SwitchTemp()
{
    return 0;
}

void GMC::set_n_temps(size_t i)
{
    if (i<2) {
        prf::cerr<<"GMC> Refusing to use less than 2 temperatures.\n";
        i=2;
    }
    make_default_temps(i);
    uph.set_n_temps(ntmp);
    explicitntmps=true;
}

void GMC::make_default_temps(size_t i)
{
    ntmp=i;
    beta.resize(i);
    dndT.resize(i,0);
    if (tmax<tmin) std::swap(tmax,tmin);
    for (unsigned c=0;c<ntmp;++c)
        beta[c]=1.0/tmax*pow(tmax/tmin,((double)c)/(ntmp-1));
}

void GMC::set_temperature_range(double t0, double t1)
{
    tmin=t0;
    tmax=t1;
    make_default_temps(ntmp);
}

int GMC::SetTemps(std::vector<double> tmptemp)
{
    if (tmptemp.size()==0) {
        prf::cerr<<"Read a total of 0 temperature values\n";
        return 0;
    }

    std::sort(tmptemp.begin(),tmptemp.end());

    ntmp=tmptemp.size();
    beta.resize(ntmp,0);
    dndT.resize(ntmp,0);

    for (unsigned i=0;i<ntmp;++i) beta[i]=tmptemp[i];
    uph.set_n_temps(ntmp);
    explicitntmps=true;
    return 1;
}

int GMC::SetTemps(std::string tempfile)
{
    if (prf_utils::TestFile_r(tempfile.c_str())==0) {
        prf::cerr<<"Reading temperatures from file"<<tempfile<<" failed.\n";
        return 0;
    }

    std::vector<std::string> lines,parts;
    std::string line;
    std::ifstream fin(tempfile.c_str());

    while (getline(fin,line)) lines.push_back(line);
    fin.close();

    bool invt=false,kelv=false;
    std::vector<double> tmparry;

    for (size_t i=0;i<lines.size();++i) {
        line=lines[i];
        line=prf_utils::trim_str(line);

        if (line.empty()) continue;

        if (line[0]=='#') {
            parts.clear();
            prf_utils::split(line,parts);

            if (parts[0]=="#temperature") {
                if ((parts[1]=="inverted") && !kelv) invt=true;
                else if (parts[1]=="Kelvin") kelv=true;
            }
        } else {
            if (prf_utils::is_number(line)) {
                double tmp=strtod(line.c_str(),NULL);
                tmparry.push_back(tmp);
            }
        }
    }

    ntmp=tmparry.size();

    if (!invt) {
        for (size_t i=0;i<tmparry.size();++i) {
            if (kelv) tmparry[i]=pru_in_kelvin/tmparry[i];
            else tmparry[i]=1.0/tmparry[i];
        }
    }
    Logger(10)<<"Read in "<<ntmp<<" temperatures from "<<tempfile<<"\n";
    if (SetTemps(tmparry)!=0) tfileinuse=true;

    return 1;
}

void GMC::writeTempStat(std::string flnm)
{
    Output ost(flnm.c_str());
    ost <<"# Number of visits to different temperatures ...\n";

    for (unsigned i=0;i<ntmp;++i) ost<<i<<"   "<<dndT[i]<<"\n";

    ost.close();
}

void GMC::writeTemperatures(std::string flnm)
{
    Output fp(flnm.c_str());
    fp<<"# Index T(model)\tT(Kelvin) \tBeta\n";

    for (unsigned i=0;i<ntmp;++i) {
        fp<<i<<"\t"<<(1.0/beta[i])<<"\t"<<(pru_in_kelvin/beta[i])
                <<"\t"<<beta[i]<<"\n";
    }

    fp.close();
}

void GMC::SetTemp(unsigned i)
{
   SetBeta(beta[itmp=i%ntmp]);
   dndT[itmp]++;
}

int GMC::prepare_continue(std::string flnm)
{
    if (prf_utils::TestFile_r(flnm)==0) {
        prf::cerr<<"Reading file "<<flnm<<" failed. Statistics on the "
                <<"number of visits to different temperatures will contain "
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

int GMC::parseCommand(InstructionString s)
{
    if (s.head()=="tfile" or s.head()=="temperature_file") {
        SetTemps(s.tail().str());
    } else if (s.head()=="ntmp" or s.head()=="num_temperatures") {
        if (!tfileinuse) set_n_temps(atoi(s.part(1).c_str()));
    } else if (s.head()=="tmax" or s.head()=="max_temperature") {
        if (!tfileinuse) {
            tmax=prf_utils::make_temperature(s.tail().str());
            set_temperature_range(tmin,tmax);
        }
    } else if (s.head()=="tmin" or s.head()=="min_temperature") {
        if (!tfileinuse) {
            tmin=prf_utils::make_temperature(s.tail().str());
            if (tmin<=0) tmin=1.0e-7;
            set_temperature_range(tmin,tmax);
        }
    }
    std::cout << "tmin="<<tmin<<"; tmax="<<tmax<<"; ntmp="<<beta.size()<<"\n";
    return MC::parseCommand(s);
}

void GMC::print_setup()
{
    prf::cout<<"\n1 Monte Carlo Cycle = "<<CycleLength()
            <<" Elementary Monte Carlo Steps or Updates\n\n";

    uph.print_setup();

    prf::cout<<"Number of Temperatures = "<<ntmp<<"\n";
}

std::string GMC::ConfSignature()
{
    std::ostringstream ost;
    ost<<"simulation_method "<<descr<<"\n";
    ost<<"ntmp "<<ntmp<<"\n";
    ost<<"temperature_array ";
    for (size_t i=0;i<beta.size();++i) {
        ost<<" "<<temperature(i);
    }
    ost<<"\n";
    return ost.str();
}
