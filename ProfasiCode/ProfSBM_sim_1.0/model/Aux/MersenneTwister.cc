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

#include "MersenneTwister.hh"
#include "prf_xml.hh"
#include <sstream>
#include <cstdlib>

using namespace prf;

static const unsigned long UPPER_MASK = 0x80000000UL;
static const unsigned long LOWER_MASK = 0x7fffffffUL;

MersenneTwister::MersenneTwister()
{
    Name("MersenneTwister");
    ResetDefaultState(0);
}

MersenneTwister::~MersenneTwister() {}

double MersenneTwister::shoot()
{
    unsigned long k ;
#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

    if (mti >= mt_N) {
        int kk;

        for (kk = 0; kk < mt_N - mt_M; kk++) {
            unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk + mt_M] ^(y >> 1) ^ MAGIC(y);
        }

        for (; kk < mt_N - 1; kk++) {
            unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk + (mt_M - mt_N)] ^(y >> 1) ^ MAGIC(y);
        }

        unsigned long y = (mt[mt_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);

        mt[mt_N - 1] = mt[mt_M - 1] ^(y >> 1) ^ MAGIC(y);
        mti = 0;
    }

    k = mt[mti];

    k ^= (k >> 11);
    k ^= (k << 7) & 0x9d2c5680UL;
    k ^= (k << 15) & 0xefc60000UL;
    k ^= (k >> 18);

    ++mti;
    ++ncalls;
    return (1.0+k) / 4294967297.0 ;
}

void MersenneTwister::ResetDefaultState(long seedd)
{
    unsigned long int s=(unsigned long int)seedd;
    int i;

    if (s == 0) s = 4357;   /* the default seed is 4357 */

    mt[0]= s & 0xffffffffUL;

    for (i = 1; i < mt_N; i++) {
        mt[i] =(1812433253UL * (mt[i-1] ^(mt[i-1] >> 30)) + i);
        mt[i] &= 0xffffffffUL;
    }

    mti = i;

    ncalls=0;
}

void MersenneTwister::saveState(std::string sttfile)
{
    Output op(sttfile.c_str());
    op<<"<?xml version=\"1.0\"?>\n<random_number_state>\n";
    op<<"<generator_type>"<<Name()<<"</generator_type>\n";
    op<<"<original_seed>"<<myseed<<"</original_seed>\n";
    op<<"<number_of_calls>"<<ncalls<<"</number_of_calls>\n";
    op<<"<data>\n";
    op<<mti<<"\n";
    for (int i=0;i<mt_N;++i) op<<mt[i]<<"\n";
    op<<"</data>\n";
    op<<"</random_number_state>\n";
    op.close();
}

int MersenneTwister::recoverState(std::string sttfile)
{
    prf_xml::XML_Node *root=prf_xml::get_xml_tree(sttfile);
     if (root==NULL) {
        prf::cerr<<Name()<<"> Could not retrieve a valid top level XML "
            <<"node from "<<sttfile<<"\n";
        return 0;
    }

    prf_xml::XML_Node * mynode=NULL;
    if (root->name()=="random_number_state") mynode=root;
    else mynode=root->child("random_number_state");
    if (mynode==NULL) {
        prf::cerr<<Name()<<"> No random_number_state found in "<<root->name()<<"\n";
        if (root) delete (root);
        return 0;
    }
    prf_xml::XML_Node *nmnode=mynode->child("generator_type");
    if (nmnode==NULL) {
        prf::cerr<<Name()<<"> Random number identification missing in "<<sttfile<<"\n";
        if (root) delete (root);
        return 0;
    }
    if (nmnode->value()!=Name()) {
        prf::cerr<<Name()<<"> The given state file was made by a generator called "
                <<nmnode->value()<<"\n"
                <<"The current generator is called "<<Name()<<"\n"
                <<"Refusing to read state information because of the name mismatch\n";
        if (root) delete (root);
        return 0;
    }
    prf_xml::XML_Node *sdnode=mynode->child("original_seed");
    if (sdnode==NULL) {
        prf::cerr<<Name()<<"> Information about the original seed is missing.\n";
        if (root) delete (root);
        return 0;
    }
    long fileseed=atol(sdnode->value().c_str());
    if (fileseed!=myseed) {
        prf::cerr<<Name()<<"> The seed found in the file was "<<fileseed<<", where as "
                <<"the present generator has been initialized with seed "<<myseed<<"\n"
                <<"Refusing to use the data in file "<<sttfile<<" because of the mismatch.\n";
        if (root) delete (root);
        return 0;
    }
    prf_xml::XML_Node *ncnode=mynode->child("number_of_calls");
    if (ncnode==NULL) {
        prf::cerr<<"Information about the number of calls is missing in file "
                <<sttfile<<". Can not determine how far the generator needs to be "
                <<"fast forwarded.\n";
        if (root) delete (root);
        return 0;
    }

    unsigned long filencalls=strtoul(ncnode->value().c_str(),NULL,10);
    if (filencalls>ncalls) {
        prf::cerr<<Name()<<"> The state stored in file "<<sttfile<<" corresponds to the time "
                <<"when the generator has made "<<filencalls<<" calls. That's in the future "
                <<"relative to the number of calls in the present generator: "<<ncalls<<"\n";
        prf::cerr<<"The data in file "<<sttfile<<" is useless to restore the state at "
                <<"ncalls.\n";
        if (root) delete (root);
        return 0;
    }
    prf_xml::XML_Node *datanode=mynode->child("data");
    if (datanode==NULL) {
        prf::cerr<<Name()<<"> Data node is empty. Failed to recover state from "<<sttfile<<"\n";
        if (root) delete (root);
        return 0;
    }
    std::istringstream ssin(datanode->value());
    ssin>>mti;
    for (int i=0;i<mt_N;++i) {
        ssin>>mt[i];
    }

    unsigned long ncallstarget=ncalls;
    ncalls=filencalls;

    prf::cout<<Name()<<"> The number of calls when the state in "<<sttfile<<" was written : "
            <<filencalls<<"\n";
    prf::cout<<Name()<<"> The number of calls recoverState() is supposed to leave behind : "
            <<ncallstarget<<"\n";
    while (ncalls<ncallstarget) shoot();
    if (ncallstarget!=filencalls) {
        prf::cout<<Name()<<"> Fast forwarded from "<<filencalls<<" to "<<ncalls
                <<" random number calls.\n";
    }
    prf::cout<<Name()<<"> Successfully restored state to the one in "<<sttfile<<"\n";

    if (root) delete (root);
    return 1;
}
