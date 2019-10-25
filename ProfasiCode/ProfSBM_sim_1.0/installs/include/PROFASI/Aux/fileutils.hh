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

#ifndef FILE_UTILS
#define FILE_UTILS

#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <valarray>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <cstdlib>


namespace prf_utils
{
    typedef enum {
        CHAR,INT,LONG,FLOAT,DOUBLE,LONG_DOUBLE, UNKNOWN
            } VarType;
    
            
    double stringtodouble(std::string st);
    //! Check if directory exists, if not create it
    void CreateDir(const char * dirname);
    //! Check if file exists and the user has permissions to read it
    int TestFile(const char * filename);
    int TestFile_r(const char * filename);
    //! Same as above, but without any error messages
    int STestFile(const char * filename);
    int STestFile_r(const char * filename);
    //! Get the contents of a file as a string
    int get_file_contents(std::string filename, std::string &buf);
    char* getFileContaining( const char* path, const char* substring, bool gIgnoreHidden );
    
    inline void CreateDir(std::string st) {CreateDir(st.c_str());}

    inline int TestFile(std::string st) {return TestFile(st.c_str());}

    inline int TestFile_r(std::string st) {return TestFile_r(st.c_str());}

    inline int STestFile(std::string st) {return STestFile(st.c_str());}

    //! Get rid of leading and trailing white space characters in a string
    /**
     * This is not really a "file" utility. But the name of this file just
     * happens to be fileutils.hh. In the future it will be probably changed
     * to misc_utils when more miscellaneous small functions accumulate here.
     */
    std::string trim_str(std::string somestr);
    //! Split a string into an array of strings based on a given character
    /**
    * Splits string into an array containing at most the first ign
    * tokens which would result because of such a spliting. The return
    * value is the number of tokens in the list of split parts.
    */
    template<class T>
    int split_str(std::string st, char spl, T &lst, unsigned int ign=1000000)
    {
        lst.clear();
        size_t i=0; unsigned ntoken=0;
        std::string tmp;

        while (i<st.size()) {
            if (st[i]==spl) {
                ++i;
                lst.push_back(tmp);tmp="";
                if (++ntoken<ign) continue; else break;
            }

            tmp.push_back(st[i++]);
        }

        if (ntoken<ign && !tmp.empty()) {lst.push_back(tmp);++ntoken;}

        return ntoken;
    }

    //! Split a string into components based on white space
    template<class T>
    void split(std::string lin, T &con)
    {
        std::string tmp="";
        size_t lgt=lin.size(),i=0;

        while (i<lgt) {
            while (isspace(lin[i])&&i<lgt) i++;

            while (!isspace(lin[i])&&i<lgt) tmp+=lin[i++];

            if (!tmp.empty()) con.push_back(tmp);

            tmp="";
        }
    }

    //! Get the contents of a file in a container of strings
    /**
      The file is read line by line, and non-empty lines are added
      to the back of a standard library container.
      */
    template<class T>
    int get_lines(std::string filename, T &lines)
    {
        if (TestFile_r(filename)) {
            lines.clear();
            std::ifstream fin(filename.c_str());
            std::string line;
            while (getline(fin,line)) {
                line=trim_str(line);
                if (!line.empty()) lines.push_back(line);
            }
            fin.close();
            return 1;
        }
        return 0;
    }
    VarType make_type(std::string);
    std::string make_type_str(int);
    bool running_on_little_endian();
    bool is_number(std::string wd0);
    double make_temperature(std::string value);
    double make_seconds(std::string value);
    //! Convert from radians to degrees
    double rad_to_deg(double r);
    //! Convert from degrees to radians
    double deg_to_rad(double d);
    std::string rm_space(std::string);
    //! Break down a filename:selection string into its components
    /**
      * Given a compound name like 1GB1.pdb::A,41,56 this function will
      * create a file name ("1GB1"), identify its file extension ("pdb") and
      * the selections (":A,41,56")
      */
    void analyze_filename(std::string inputnm, std::string &filename, std::string &extn, std::string &sel);
    void make_rel_path(std::string flnm, std::list<std::string> &pth);
}

#endif
