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

#include "RandomNumberBase.hh"
#include "fileutils.hh"
#include <sstream>
#include <ctime>
#include <sys/time.h>
#include "profasi_io.hh"
#include <sys/types.h>
#include <unistd.h>

using namespace prf;
RandomNumberBase::RandomNumberBase() :myseed(-11)
{
    Name("Random_Number_Base");
    ncalls=0;
}

RandomNumberBase::~RandomNumberBase() {}

double RandomNumberBase::shoot() {++ncalls; return 0.5;}

void RandomNumberBase::WriteTo(FILE *cnfil)
{
    fwrite(&myseed,sizeof(long),1,cnfil);
    fwrite(&ncalls,sizeof(unsigned long),1,cnfil);
}

void RandomNumberBase::ReadFrom(FILE *cnfil)
{
    fread(&myseed,sizeof(long),1,cnfil);
    fread(&ncalls,sizeof(unsigned long),1,cnfil);
}

void RandomNumberBase::WriteTo_text(FILE *cnfil)
{
    fprintf(cnfil,"%li\n%lu\n",myseed,ncalls);
}

void RandomNumberBase::ReadFrom_text(FILE *cnfil)
{
    fscanf(cnfil,"%li\n%lu\n",&myseed,&ncalls);
}

void RandomNumberBase::ResetDefaultState(long seedd) {ncalls=0;}

void RandomNumberBase::RandomizeState(int inpseed)
{
    myseed=createSeed(inpseed);
    ResetDefaultState(myseed);
}

int RandomNumberBase::createSeed(int inpseed)
{
    int i1=getpid();
    timeval tnow;
    gettimeofday(&tnow,NULL);

    return tnow.tv_usec-i1*(inpseed+1);
}

void RandomNumberBase::setSeed(long i)
{
    myseed=i;
    ResetDefaultState(myseed);
}

int RandomNumberBase::ConfSize()
{
    return (sizeof(long)+sizeof(unsigned long))/sizeof(char);
}

std::string RandomNumberBase::ConfSignature()
{
    std::ostringstream sout;
    sout<<"long 1 "<<sizeof(long)/sizeof(char)
    <<" "<<sizeof(long)/sizeof(char)<<"\n"
    <<"unsigned long 1 "<<sizeof(unsigned long)/sizeof(char)
    <<" "<<sizeof(unsigned long)/sizeof(char)<<"\n";
    return sout.str();
}

void RandomNumberBase::retrace(std::string sttfile)
{
    Logger blog(5);
    if (sttfile.empty() || (prf_utils::TestFile_r(sttfile)==0)
        || recoverState(sttfile)==0) {
        prf::cout<<Name()<<"> Resetting seed to "<<myseed<<"\n";
        unsigned long ncalls0=ncalls;
        ResetDefaultState(myseed);
        prf::cout<<Name()<<"> Retracing generation history through "<<ncalls0<<" calls.\n";

        for (unsigned long i=0;i<ncalls0;++i) shoot();

        prf::cout<<Name()<<"> Done!\n";
    }
}

void RandomNumberBase::saveState(std::string sttfile)
{
    Output op(sttfile.c_str());
    op<<"<?xml version=\"1.0\"?>\n<random_number_state>\n";
    op<<"<generator_type>"<<Name()<<"</generator_type>\n";
    op<<"<original_seed>"<<myseed<<"</original_seed>\n";
    op<<"<number_of_calls>"<<ncalls<<"</number_of_calls>\n";
    op<<"<data></data>\n";
    op<<"</random_number_state>\n";
    op.close();
}

int RandomNumberBase::recoverState(std::string sttfile)
{
    prf::cerr<<"Can not recover state from "<<sttfile<<" using base class "
            <<"recover function! Seems like recoverState() has not been "
            <<"implemented for "<<Name()<<"!\n";
    return 0;
}
