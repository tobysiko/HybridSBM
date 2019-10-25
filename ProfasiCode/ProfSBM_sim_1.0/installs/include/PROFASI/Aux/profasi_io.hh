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

#ifndef profasi_io_HH
#define profasi_io_HH
#include <cstdio>
#include <string>
#include "Vector3.hh"

namespace prf_utils
{

    class Seconds
    {
    public:
        explicit Seconds(double);
        ~Seconds();
        inline double val() const {return v;}

    private:
        double v;
    };
}

namespace prf
{

    class Output
    {
    private:
        FILE *mh;
    public:
        Output();
        Output(const char *filename, const char * mode="w");
        Output(FILE *filehandl);
        Output(const Output &);
        ~Output();
        Output &operator=(const Output &);

        void Attach(const char *filename, const char * mode="w");
        void Attach(FILE *filehandl);
        void open(const char *filename, const char * mode="w");
        void close();
        void flush();
        void flush(std::string s);
        Output & operator<<(int i);
        Output & operator<<(long i);
        Output & operator<<(double x);
        Output & operator<<(unsigned int i);
        Output & operator<<(unsigned long i);
        Output & operator<<(float i);
        Output & operator<<(char ch);
        Output & operator<<(const char *ch);
        inline Output & operator<<(std::string s) {return (operator<<(s.c_str()));}

        Output & operator<<(const Vector3 &v);
        Output & operator<<(const prf_utils::Seconds &s);
    };

    class Logger
    {
    public:
        Logger();
        Logger(Output *);
        Logger(int ilevel);
        ~Logger();
        Logger(const Logger &);
        Logger & operator=(const Logger &);
        void flush();
        Logger & operator()(int msg_level);
        inline int current_status() const {return current_level;}

        inline void verbosity_level(int thrs) {verbosity=thrs;}

        template<typename T>
        Logger & operator<<(T x) {
            if (verbosity>=current_level)(*op)<<x;

            return *this;
        }

    private:
        int current_level;
        Output *op;
    public:
        static int verbosity;
    };

    void greet();
    void highlight(std::string sr="");
    inline void highlight(const char *cst){highlight(std::string(cst));}

    extern Output clog;
    extern Output cout;
    extern Output cerr;
}

#endif
