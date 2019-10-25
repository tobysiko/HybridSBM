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

#include "SimAnneal.hh"

SimAnneal::SimAnneal()
{
    schemafile="none.sas";
    nswatT=500;
    useschemafile=false;
    par.option("schema_file", "sch", 1,"(simulated annealing schema)");
    par.option("cycles_at_T","ncT",1,
                "(in absence of explicit schema file)");
    descr="Simulated annealing";
}

SimAnneal::~SimAnneal() {}

int SimAnneal::Setup()
{
    if (useschemafile) useschemafile=(read_schema()==1);
    int nerr=MC::Setup();
    if (!useschemafile) make_default_schema();

    counter=itmp=0;
    SetTemp(itmp);
    return nerr;
}

unsigned SimAnneal::SwitchTemp()
{
    if (++counter==ncbeta[itmp]) {
        ++itmp;
        counter=0;
    }

    if (itmp==ntmp) counter=0; else SetTemp(itmp);

    return itmp;
}

void SimAnneal::print_setup()
{
    GMC::print_setup();
    prf::cout<<"Index\tbeta\t\tT.(model)\tT.(kelvin)\t cycles at T \n";

    for (size_t i=0;i<ntmp;++i) {
        prf::cout<<i<<'\t'<<beta[i]<<'\t'<<(1.0/beta[i])<<'\t'
                <<(UnivConstants::pru_in_kelvin/beta[i])<<'\t'<<ncbeta[i]<<'\n';
    }
}

int SimAnneal::parseCommand(InstructionString s)
{
    if (s.head()=="ncT" or s.head()=="cycles_at_T") {
        nswatT=atoi(s.tail().str().c_str());

        if (nswatT<=0) {
            prf::cerr<<"Inappropriate value "<<nswatT
                    <<" given for sweeps at each temperature\n";
            prf::cerr<<"Using default, 500 sweeps instead\n";
            nswatT=500;
        }
    } else if (s.head()=="schema_file") {
        schemafile=s.tail().str();
        useschemafile=true;
    }
    return GMC::parseCommand(s);
}


int SimAnneal::read_schema()
{
    if (prf_utils::TestFile_r(schemafile.c_str())==0) return 0;

    std::vector<std::string> lines,parts;
    std::string line;
    std::ifstream fin(schemafile.c_str());

    while (getline(fin,line)) lines.push_back(line);
    fin.close();

    bool invt=false,kelv=false;
    beta.clear();
    ncbeta.clear();

    for (size_t i=0;i<lines.size();++i) {
        parts.clear();
        line=lines[i];
        line=prf_utils::trim_str(line);

        if (line.empty()) continue;

        prf_utils::split(line,parts);

        if (line[0]=='#') {
            if (parts[0]=="#temperature") {
                if (parts[1]=="inverted" && !kelv) invt=true;
                else if (parts[1]=="Kelvin") kelv=true;
            }
        } else {
            if (prf_utils::is_number(parts[0]) &&
                prf_utils::is_number(parts[1])) {
                double tmp=strtod(parts[0].c_str(),NULL);
                int frq=atoi(parts[1].c_str());
                beta.push_back(tmp);
                ncbeta.push_back(frq);
            }
        }
    }

    if (beta.size()<2) return 0; else {
        ntmp=beta.size();
        uph.set_n_temps(ntmp);
        dndT.resize(ntmp,0);
        explicitntmps=true;
    }
    if (not invt) {
        for (size_t i=0;i<beta.size();++i) {
            if (kelv) beta[i]=UnivConstants::pru_in_kelvin/beta[i];
            else beta[i]=1.0/beta[i];
        }
    }
    Logger(10)<<"Read simulated annealing schema from file "<<schemafile<<"\n";
    return 1;
}

int SimAnneal::make_default_schema()
{
    beta.resize(ntmp,0);
    ncbeta.resize(ntmp,nswatT);
    double tfact=pow((tmax/tmin),(1.0/(ntmp-1)));

    for (size_t i=0;i<ntmp;++i) beta[i]=(1.0/tmax)*pow(tfact,(double) i);

    return 1;
}

