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

#include "profasi_io.hh"
#include "profasi_version.hh"
#include <time.h>
#include <unistd.h>
#include "prf_time.hh"
using std::string;

using namespace prf_utils;

namespace prf
{

    Output::Output() : mh(NULL) {}

    Output::Output(const char *filename,const char * mode) : mh(NULL)
    {
        mh=fopen(filename,mode);
    }

    Output::Output(FILE *filehandl) : mh(NULL)
    {
        mh=filehandl;
    }

    Output::~Output() {close();}

    Output::Output(const Output &op) : mh(NULL)
    {
        mh=op.mh;
    }

    Output &Output::operator=(const Output &op)
    {
        if (this!=&op) {
            mh=op.mh;
        }

        return *this;
    }

    void Output::Attach(const char *filename,const char * mode)
    {
        mh=fopen(filename,mode);
    }

    void Output::Attach(FILE *filehandl)
    {
        mh=filehandl;
    }

    Output & Output::operator<<(int i)
    {
        fprintf(mh,"%d",i);
        return *this;
    }

    Output & Output::operator<<(long i)
    {
        fprintf(mh,"%li",i);
        return *this;
    }

    Output & Output::operator<<(double x)
    {
        fprintf(mh,"%f",x);
        return *this;
    }

    Output & Output::operator<<(unsigned int i)
    {
        fprintf(mh,"%d",i);
        return *this;
    }

    Output & Output::operator<<(unsigned long i)
    {
        fprintf(mh,"%lu",i);
        return *this;
    }

    Output & Output::operator<<(float x)
    {
        fprintf(mh,"%f",x);
        return *this;
    }

    Output & Output::operator<<(char ch)
    {
        fprintf(mh,"%c",ch);
        return *this;
    }

    Output & Output::operator<<(const char *ch)
    {
        fprintf(mh,"%s",ch);
        return *this;
    }

    Output & Output::operator<<(const Vector3 &v)
    {

        fprintf(mh,"%f\t%f\t%f",v.x(),v.y(),v.z());
        return *this;
    }

    Output &Output::operator<<(const Seconds &s)
    {
        int days=0,hours=0,minutes=0;
        double v=s.val();
        string tunits=" seconds";

        if (v>60) {
            minutes=(int)(v/60);
            v-=60*minutes;
        }

        if (minutes>=60) {
            hours=minutes/60;
            minutes=minutes%60;
        }

        if (hours>=24) {
            days=hours/24;
            hours=hours%24;
        }

        if (days>0) {
            operator<<(days)<<" days ";
        }

        if (hours>0) {
            operator<<(hours)<<" hours ";
        }

        if (minutes>0) {
            operator<<(minutes)<<" minutes "
            <<((int)v)<<" seconds ";
        } else fprintf(mh,"%4.2f seconds ",v);

        return *this;
    }

    void Output::open(const char * filename,const char * mode)
    {
        if (mh!=NULL && mh!=stdout && mh!=stderr) fclose(mh);

        mh=fopen(filename,mode);
    }

    void Output::close()
    {
        if (mh!=NULL && mh!=stdout && mh!=stderr) {
            fclose(mh);
            mh=NULL;
        }
    }

    void Output::flush()
    {
        fflush(mh);
    }

    void Output::flush(std::string s)
    {
        fprintf(mh,"%s",s.c_str());
        fflush(mh);
    }

    Logger::Logger()
    {
        op=&prf::clog;
        current_level=1;
    }

    Logger::Logger(Output *opt)
    {
        op=opt;
        current_level=1;
    }

    Logger::Logger(int ilevel)
    {
        op=&prf::clog;
        current_level=ilevel;
    }

    Logger::~Logger() {}

    Logger::Logger(const Logger &lg)
    {
        op=lg.op;
        current_level=lg.current_level;
    }

    Logger & Logger::operator=(const Logger &lg)
    {
        if (this!=&lg) {
            op=lg.op;
            current_level=lg.current_level;
        }

        return *this;
    }
    void Logger::flush()
    {
        op->flush();
    }

    Logger &Logger::operator()(int new_level)
    {
        current_level=new_level;
        return *this;
    }

    void greet()
    {
        prf::cout<<"\n\n\n";
        highlight(("ProFASi v."+profasi_version()).c_str());
        prf::cout<<"ProFASi:  Protein Folding and Aggregation Simulator.\n"
        <<"Copyright (2005) Anders Irback and Sandipan Mohanty\n"
        <<"Reference: J. Comput. Chem. [27] 1548 (2006)\n"
        <<"Web address: http://cbbp.thep.lu.se/activities/profasi/\n";
        prf::cout<<"Version details "<<profasi_git_version()<<"\n";
        char hostname[30];
        prf::cout<<"Starting at UTC "<<prf_time().to_UTC()<<"\n";
        prf::cout.flush();
        gethostname(hostname,30);
        prf::cout<<"Hosting machine is "<<hostname<<"\n";
        prf::cout.flush();
        highlight();
        prf::cout.flush();
    }

    void highlight(string blaha)
    {
        int fl=30-blaha.size()/2;

        for (int i=0;i<fl;++i) prf::cout<<'>';prf::cout<<' ';

        prf::cout<<blaha;prf::cout<<' ';

        for (int i=0;i<fl;++i) prf::cout<<'<';prf::cout<<'\n';
    }

    Output cout(stdout);
    Output cerr(stderr);
    Output clog(stdout);
    int Logger::verbosity=100;
}

namespace prf_utils
{
    Seconds::Seconds(double tx) {v=tx;}

    Seconds::~Seconds() {}
}
