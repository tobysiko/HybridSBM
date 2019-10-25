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

#include "fileutils.hh"
#include "profasi_io.hh"
#include "Constants.hh"
#include <vector>
#include <list>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <dirent.h>
//#include <boost/regex.hpp>
//#include <boost/algorithm/string/replace.hpp>

using namespace UnivConstants;
using std::string;

using namespace prf;

namespace prf_utils
{


/**
	bool MatchTextWithWildcards(const string &text, string wildcardPattern, bool caseSensitive )//= true)
	{
	    // Escape all regex special chars
	    EscapeRegex(wildcardPattern);
	
	    // Convert chars '*?' back to their regex equivalents
	    boost::replace_all(wildcardPattern, "\\?", ".");
	    boost::replace_all(wildcardPattern, "\\*", ".*");
	
	    boost::regex pattern(wildcardPattern, caseSensitive ? regex::normal : regex::icase);
	
	    return regex_match(text, pattern);
	}
	
	void EscapeRegex(string &regex)
	{
	    boost::replace_all(regex, "\\", "\\\\");
	    boost::replace_all(regex, "^", "\\^");
	    boost::replace_all(regex, ".", "\\.");
	    boost::replace_all(regex, "$", "\\$");
	    boost::replace_all(regex, "|", "\\|");
	    boost::replace_all(regex, "(", "\\(");
	    boost::replace_all(regex, ")", "\\)");
	    boost::replace_all(regex, "[", "\\[");
	    boost::replace_all(regex, "]", "\\]");
	    boost::replace_all(regex, "*", "\\*");
	    boost::replace_all(regex, "+", "\\+");
	    boost::replace_all(regex, "?", "\\?");
	    boost::replace_all(regex, "/", "\\/");
	}*/
	char* getFileContaining( const char* path, const char* substring, bool gIgnoreHidden )
	{
	
	   DIR* dirFile = opendir( path );
	   if ( dirFile ) 
	   {
	      struct dirent* hFile;
	      //errno = 0;
	      while (( hFile = readdir( dirFile )) != NULL ) 
	      {
	         if ( !strcmp( hFile->d_name, "."  )) continue;
	         if ( !strcmp( hFile->d_name, ".." )) continue;

	         // in linux hidden files all start with '.'
	         if ( gIgnoreHidden && ( hFile->d_name[0] == '.' )) continue;

	         // dirFile.name is the name of the file. Do whatever string comparison 
	         // you want here. Something like:
	         if ( strstr( hFile->d_name, substring )){
	        	 printf( "found an %s file: %s", substring, hFile->d_name );
	        	 closedir( dirFile );
	        	 return hFile->d_name;
	         }
	            
	         	
	      } 
	      closedir( dirFile );
	   }
	   return "null";
	}
	
	double stringtodouble (std::string val){
		//return std::stod(val);
		double t;
		std::stringstream convert( val  );
		convert>>t;
		return t;
	}
    void CreateDir(const char * dirname)
    {
        mkdir(dirname,(mode_t) 0755);
    }

    int TestFile(const char * filename)
    {
        if (access(filename,F_OK)==-1) {
            prf::cerr<<"File \""<<filename<<"\" not found. \n";
            return 0;
        }

        if (access(filename, R_OK|W_OK)!=0) {
            prf::cerr<<"access denied to file \""<<filename<<"\".\n";
            return 0;
        }

        return 1;
    }

    int TestFile_r(const char * filename)
    {
        if (access(filename,F_OK)==-1) {
            prf::cerr<<"File \""<<filename<<"\" not found. \n";
            return 0;
        }

        if (access(filename, R_OK)!=0) {
            prf::cerr<<"access denied to file \""<<filename<<"\".\n";
            return 0;
        }

        return 1;
    }

    int STestFile(const char * filename)
    {
        if (access(filename,F_OK)==-1) {
            return 0;
        }

        if (access(filename, R_OK|W_OK)!=0) {
            return 0;
        }

        return 1;
    }

    int STestFile_r(const char * filename)
    {
        if (access(filename,F_OK)==-1) {
            return 0;
        }

        if (access(filename, R_OK)!=0) {
            return 0;
        }

        return 1;
    }

    int get_file_contents(std::string flname, std::string &buf)
    {
        if (TestFile_r(flname.c_str())==0) return 0;
        FILE *fp=fopen(flname.c_str(),"r");
        char c;
        while ((c=getc(fp))!=EOF) buf.push_back(c);
        fclose(fp);
        return buf.size();
    }

    string trim_str(string somestr)
    {
        if (somestr.empty()) return somestr;
        string trmmd="";
        size_t i=0,j=somestr.size()-1;

        while (isspace(somestr[i]) && i<=j) ++i;

        while (isspace(somestr[j]) &&j>=i) --j;

        return somestr.substr(i,j-i+1);
    }

    bool running_on_little_endian()
    {
        //Ref: www.faqs.org/faqs/graphics/fileformats-faq/part4/section-7.html
        short int word = 0x0001;
        char *byte = (char *) &word;
        return(byte[0] ? true : false);
    }

    bool is_number(string wd0)
    {
        string wd;
        bool ans=true;
        int ndeci=0,strt=0,iend=(int) wd0.size();

        if (!iend) return false;

        while (isspace(wd0[strt])) ++strt;

        while (isspace(wd0[iend-1])) --iend;

        wd=wd0.substr(strt,iend-strt);

        for (size_t i=0;ans&&i<wd.size();++i) {
            if (wd[i]=='-' && i!=0) ans=false;

            if (wd[i]=='.') ++ndeci;

            if (!(isdigit(wd[i])||wd[i]=='-'||wd[i]=='.')) ans=false;
        }

        if (ndeci>1) ans=false;

        return ans;
    }

    double make_temperature(string value)
    {
        std::vector<string> parts;
        split<std::vector<string> >(value,parts);
        double tvl=strtod(parts[0].c_str(),NULL);
        double convfact=1.0;

        if (parts.size()>1) {
            if (parts[1]==string("K") || parts[1]==string("Kelvin")) {
                convfact=kelvin_in_pru;
            }
        }

        return convfact*tvl;

    }

    double deg_to_rad(double d)
    {
        const double convfact=UnivConstants::pi/180.0;
        return convfact*d;
    }

    double rad_to_deg(double r)
    {
        const double convfact=180/UnivConstants::pi;
        return convfact*r;
    }

    double make_seconds(std::string value)
    {
        std::list<string> parts;
        split_str<std::list<string> >(value,':',parts,4);
        const double convs[]={1,60,3600,86400};
        double ans=0;
        int j=0;

        for (std::list<string>::reverse_iterator i=parts.rbegin();
             i!=parts.rend();++i,++j) {
            ans+=convs[j]*strtod(i->c_str(),NULL);
        }

        return ans;
    }

    std::string rm_space(std::string inp)
    {
        size_t i0=0,i1=inp.size();
        while (i0<i1 && isspace(inp[i0])) ++i0;
        while (i1>i0 && isspace(inp[i1-1])) --i1;
        return std::string(inp,i0,i1-i0);
    }

    void analyze_filename(std::string ifile, std::string &filename, std::string &extn, std::string &sel)
    {
        size_t icolon=ifile.find(':');
        if (icolon<ifile.size()-1) sel=std::string(ifile,icolon+1);
        else sel="1:*";
        filename=std::string(ifile,0,icolon);
        size_t idot=filename.find_last_of('.');
        extn=string(filename,idot+1);
    }

    VarType make_type(std::string tstr)
    {
        if (tstr=="int") return INT;
        if (tstr=="long") return LONG;
        if (tstr=="char") return CHAR;
        if (tstr=="float") return FLOAT;
        if (tstr=="double") return DOUBLE;
        if (tstr=="long double") return LONG_DOUBLE;
        if (tstr=="unsigned int") return INT;
        if (tstr=="unsigned long") return LONG;
        if (tstr=="unsigned char") return CHAR;
        return UNKNOWN;
    }

    std::string make_type_str(int vt)
    {
        std::string ans="UNKNOWN";
        switch(vt) {
        case INT: ans="int"; break;
        case LONG: ans="long"; break;
        case CHAR: ans="char"; break;
        case FLOAT: ans="float";break;
        case DOUBLE: ans="double"; break;
        case LONG_DOUBLE: ans="long double"; break;
        default: ans="UNKNOWN"; break;
        };
        return ans;
    }

    void make_rel_path(std::string flnm, std::list<std::string> &pth)
    {
        pth.clear();
        prf_utils::split_str(flnm,'/',pth);
        std::list<std::string>::iterator it=pth.begin(),jt=pth.begin();
        while  (jt!=pth.end()) {
            ++jt;
            std::string d=*it;
            if (d.empty() or d==".") pth.erase(it);
            else if (d==".." and it!=pth.begin()) {
                std::list<std::string>::iterator kt=it;
                --kt;
                if (*kt!="..") {
                    pth.erase(kt);
                    pth.erase(it);
                }
            }
            it=jt;
        }
    }
}
