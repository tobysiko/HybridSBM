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

#include <Aux/ProgUtils.hh>
#include <Aux/fileutils.hh>
#include <Aux/profasi_io.hh>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <zlib.h>

#define bufsize 1024

using namespace prf_utils;
using std::vector;
using std::string;
using std::ifstream;
using std::istringstream;
using std::deque;

void typeseq_append(int n, string tstr, vector<string> & tslst)
{
    for (int j=0;j<n;++j) tslst.push_back(tstr);
}

bool splice_type_string(deque<string> &lparts, string &typstr)
{
    bool ans=false;
    vector<string> knowntypes;
    knowntypes.push_back("int");
    knowntypes.push_back("char");
    knowntypes.push_back("long");
    knowntypes.push_back("float");
    knowntypes.push_back("double");

    typstr="";

    if (lparts.empty()) return false;

    if (lparts.front()==string("unsigned")) {
        typstr="unsigned ";
        lparts.pop_front();
    }

    if (find(knowntypes.begin(),
             knowntypes.end(),lparts.front())!=knowntypes.end()) {
        typstr+=lparts.front();
        ans=true;
    }

    string lstkwd=lparts.front();

    if (ans) lparts.pop_front();

    if (lstkwd==string("long")) {
        if (find(knowntypes.begin(),
                 knowntypes.end(),(*lparts.begin()))!=knowntypes.end()) {
            typstr+=(*lparts.begin());
            lparts.pop_front();
        }
    }

    return ans;
}

bool get_next_var(FILE *fp, string typstr, string &res)
{
    char ans[20];
    size_t len=0;

    if (typstr==string("int")) {
        int i;
        len=fread(&i,sizeof(int),1,fp);
        sprintf(ans,"%d",i);
    } else if (typstr==string("long")) {
        long i;
        len=fread(&i,sizeof(long),1,fp);
        sprintf(ans,"%li",i);
    } else if (typstr==string("char")) {
        char i;
        len=fread(&i,sizeof(char),1,fp);
        sprintf(ans,"%c",i);
    } else if (typstr==string("float")) {
        float i;
        len=fread(&i,sizeof(float),1,fp);
        sprintf(ans,"%f",i);
    } else if (typstr==string("double")) {
        double i;
        len=fread(&i,sizeof(double),1,fp);
        sprintf(ans,"%.16f",i);
    } else if (typstr==string("unsigned int")) {
        unsigned int i;
        len=fread(&i,sizeof(unsigned int),1,fp);
        sprintf(ans,"%u",i);
    } else if (typstr==string("unsigned long")) {
        unsigned long i;
        len=fread(&i,sizeof(unsigned long),1,fp);
        sprintf(ans,"%lu",i);
    } else if (typstr==string("unsigned char")) {
        unsigned char i;
        len=fread(&i,sizeof(unsigned char),1,fp);
        sprintf(ans,"%uc",i);
    }

    res=string(ans);

    return (len!=0);
}

bool put_next_var(string locvar, string typstr, FILE *fp)
{
    size_t len=0;

    if (typstr==string("int")) {
        int i=atoi(locvar.c_str());
        len=fwrite(&i,sizeof(int),1,fp);
    } else if (typstr==string("long")) {
        long i=atol(locvar.c_str());
        len=fwrite(&i,sizeof(long),1,fp);
    } else if (typstr==string("char")) {
        char i=locvar[0];
        len=fwrite(&i,sizeof(char),1,fp);
    } else if (typstr==string("float")) {
        float i=atof(locvar.c_str());
        len=fwrite(&i,sizeof(float),1,fp);
    } else if (typstr==string("double")) {
        double i=strtod(locvar.c_str(),NULL);
        len=fwrite(&i,sizeof(double),1,fp);
    } else if (typstr==string("unsigned int")) {
        unsigned int i=(unsigned) atoi(locvar.c_str());
        len=fwrite(&i,sizeof(unsigned int),1,fp);
    } else if (typstr==string("unsigned long")) {
        unsigned long i=strtoul(locvar.c_str(),NULL,10);
        len=fwrite(&i,sizeof(unsigned long),1,fp);
    } else if (typstr==string("unsigned char")) {
        unsigned char i=(unsigned char) locvar[0];
        len=fwrite(&i,sizeof(unsigned char),1,fp);
    }

    return (len!=0);
}

bool read_block(FILE *fp, vector<string> &tseq, string &confblock)
{
    string dummy;
    size_t ivr;
    confblock.clear();

    for (ivr=0;ivr<tseq.size();++ivr) {
        if (get_next_var(fp,tseq[ivr],dummy)) {
//            prf::cout<<"next var = "<<dummy<<"\n";
            confblock+=(dummy+"\n");
        } else break;
    }

    return (ivr==tseq.size());
}

void bufcopy(string blk, char *gzipbuf, size_t &gbufloc, gzFile gp)
{
    size_t i=0;

    do {
        while (gbufloc<bufsize && i< blk.size()) {
            gzipbuf[gbufloc++]=blk[i++];
        }

//        prf::cout<<"feeding stopped at i="<<i<<" and gbufloc = "<<gbufloc<<"\n";
        if (gbufloc==bufsize) {
            gzwrite(gp,gzipbuf,(unsigned)bufsize);
//            prf::cout<<"wrote to zip file and reset location\n";
            gbufloc=0;
        }
    } while (i<blk.size());
}

void flushgzbuf(char *gzbuf, size_t bufloc, gzFile gp)
{
//    prf::cout<<"flushgzbuf: writing the remaining "<<bufloc<<" bytes\n";
//    prf::cout<<"the current state of the gzip buffer is ...\n";
    //for (int i=0;i<bufsize;++i) printf("%c",gzbuf[i]);
    gzwrite(gp, gzbuf, bufloc);
}

bool get_next_gblock(gzFile in, string &newbuf)
{
    char gzbuf[bufsize];
    newbuf.clear();
    int len = gzread(in, (char *) gzbuf, (sizeof(gzbuf)));
//    prf::cout<<"get_next_gblock: read in "<<len<<" bytes\n";

    if (len<=0) return false;

    newbuf=string(gzbuf,len);

    return true;
}

int main(int argc, char *argv[])
{
    ProgArgs par;
    par.option("output_file","o",1);
    par.option("info_file","i",1);
    par.option("reverse","r",0,"(make a conf file from gzipped text conf)");
    par.analyze(argc,argv);

    if (par.n_spare_args()==0) {
        prf::cout<<"Usage: \n\n";
        prf::cout<<argv[0]<<" [parameters] conf_file\n\n";
        prf::cout<<"parameters could be any number and combination of ...\n";
        par.write_available();
        prf::cout<<"Example: \n\n";
        prf::cout<<argv[0]<<" -i n0/conf.info -o gztconf n0/conf\n\n";
        return 0;
    }

    string ifile=par.spare_args(0), cinfo="conf.info",line,typ,wd;

    if (par.option_given("i")) cinfo=par.option("i");
    else if (STestFile_r(cinfo.c_str())==0) cinfo=ifile+".info";

    if (TestFile_r(ifile.c_str())==0 || TestFile_r(cinfo.c_str())==0) return 1;

    string ofile="gztconf";

    if (par.option_given("o")) ofile=par.option("o");

    ifstream fin(cinfo.c_str());

    vector<string> lines, tseq;

    while (getline(fin,line)) lines.push_back(line);

    int cblocksize=bufsize;

    for (size_t i=0;i<lines.size();++i) {
        line=lines[i];
        deque<string> parts;
        split<deque<string> >(line,parts);
//        if (parts.front()==string("conf_length"))
//            cblocksize=atoi(parts[1].c_str());
        int noftyp=0;

        if (splice_type_string(parts,typ)) {
            noftyp=atoi((*parts.begin()).c_str());
            typeseq_append(noftyp,typ,tseq);
        }
    }

//    prf::cout<<"Determined conf file layout as ...\n";
//    for (size_t i=0;i<tseq.size();++i) {
//        prf::cout<<tseq[i]<<"\n";
//    }
    if (!par.option_given("r")) {
        char *gzbuf=new char[cblocksize];
        FILE *fp=fopen(ifile.c_str(),"r");
        gzFile gp=gzopen(ofile.c_str(),"w");
        string blockbuf;
        size_t gzbufloc=0;

        while (read_block(fp,tseq,blockbuf)) {
//            prf::cout<<"Read conf block\n";
//            prf::cout<<"character block size = "<<blockbuf.size()<<"\n";
            bufcopy(blockbuf,(char *)gzbuf,gzbufloc, gp);
        }

        flushgzbuf((char *)gzbuf,gzbufloc,gp);

        fclose(fp);
        gzclose(gp);
        delete(gzbuf);
        return 0;
    } else {
        char *gzbuf=new char[cblocksize];
        gzFile gp=gzopen(ifile.c_str(),"r");
        FILE *fp=fopen(ofile.c_str(),"w");
        string textconf,tmpstr,nxt;
        size_t ti=0,ncomplete=0;
        bool input_open=get_next_gblock(gp,textconf);
        istringstream *ssin= new istringstream(textconf);
        size_t remchars=textconf.size();

        while ((*ssin)>>tmpstr) {
            remchars-=(tmpstr.size()+1);
//            prf::cout<<"working on string "<<tmpstr<<" remchars = "
//                    <<remchars<<"\n";
            put_next_var(tmpstr, tseq[ti++], fp);

            if (ti==tseq.size()) {
                ++ncomplete;
                ti=0;
            }

            if (remchars<100 && input_open) {
                input_open=get_next_gblock(gp,nxt);
                textconf=textconf.substr(textconf.size()-remchars)+nxt;
                delete(ssin);
                ssin=new istringstream(textconf);
                remchars=textconf.size();
            }
        }

        fclose(fp);

        gzclose(gp);
        delete(gzbuf);

        if (ssin) delete(ssin);

        return 0;
    }
}

/**
\page conf_to_gzip Converting binary conf files to gzipped text and back
The program conf_to_gzip can convert back and forth between the binary
 configuration file layout from PROFASI simulation programs and a gzipped text
 format.

To continue a run performed on machine A on machine B with a different
 architecture, or to analyze data generated in machine A on machine B, it is
 often necessary to browse the simulation history as stored in the configuration
 files. But this poses an immediate problem: the binary files generated on one
 machine, when read on another, give rise to nonsense. This program is meant to
 deal with precisely this situation.

There are three steps involved:
<ol>
<li>Use conf_to_gzip to convert the binary format to a compressed gzipped
 format of the same data. Of course, we are not gzipping our binary data. That
 would not really "compress" anything. The binary data is already compressed
 taking advantage of apriori knowledge of the layout of the data. This can not be
 further compressed with a general purpose algorithm like the "deflate" algorithm
 of gzip. conf_to_gzip.ex interprets the binary data as they would be interpreted
 in a simulation, and then creates a gzipped file of the text representation of that data. So, if you gunzip the file generated by the program, you will get
 an ordinary text file.  </li>
<li>Copy the gzip files to your target system. As described above, these are
 gzipped versions of text files. So, you can transfer them anywhere, and they
 will still mean the same thing. </li>
<li>Use conf_to_gzip on the target machine to convert the gzip file to a
 binary configuration file suitable for that machine. </li>
</ol>
To interpret the binary format, the program needs a map of the data layout in
 that particular conf file. This information is generated by the same simulation
 programs, and stored in the same directory where they store their "conf" files.
 These files are named "conf.info". The program conf_to_gzip requires the
 location of a "conf.info" file to be able to interpret the conf data correctly.

\section confconvert Usage example

Basic usage example:\n\n
<tt>
$ conf_to_gzip -i conf.info conf -o conf.txt.gz
</tt>

This takes an input binary configuration file called conf, a file called "conf.info" containing a map of the binary data in conf and produces an output file in gzipped format called "conf.txt.gz". The output file is a gzipped text file. On most modern systems, one can simply "less" it and view its contents without uncompressing it. For the reverse process, making a binary file out of a gzipped text configuration file,

<tt>
$ conf_to_gzip -r -i conf.info conf.txt.gz -o conf
</tt>

Note that the option "i" does not stand for input. It is "info_file". The input in this second case is conf.txt.gz, and output is the binary file conf.

You have run ParTemp.mex on a cluster for 3 weeks. While analyzing the data on
 your laptop, you realize that node 13 had visited a very interesting state on
 cycle 1599999, but that state was not the lowest energy state. You want to use
 extract_snapshot to get a PDB of that state and see what it looks like. You
 need the program state histories saved in files n0/conf, n1/conf etc for each
 compute node. But the cluster consists of AMD Opteron processors, while your
 laptop has a 32 bit Intel Pentium M processor. The binary conf files are
 incompatible.


Of course, one solution is to log in to the cluster, and run the
 extract_snapshot command there. But if this is not possible for some reason,
 you should proceed as follows...

Whenever you copy your run data from one system to another, you should remember
 to convert the "conf" files to the native format of the new system. So, before
 you copy anything, run the following command in the directory containing n0, n1
 etc. ...

for j in `seq 0 max_run_index` ; do <br>
conf_to_gzip -i n$j/conf.info n$j/conf -o n$j/conf.txt.gz <br>
done<br>

This should create gzipped text configurations named "conf.txt.gz" in each
 directory. You might want to add another line "rm -f n$j/conf" before the
"done" statement to remove the binary conf file. It can be regenerated when
 required from the conf.txt.gz file when required on any system.

Next, copy your run directory to your laptop. Now you can to create binary conf
 files for the laptop, as if they had been natively generated there. Type the
 following:

for j in `seq 0 max_run_index`; do <br>
conf_to_gzip -r -i n$j/conf.info n$j/conf.txt.gz -o n$j/conf<br>
done<br>

Now, you can proceed with the structure extraction with "extract_snapshot" in
 the usual way. This same procedure can be applied to continue one set of runs
 done on one cluster on another.

*/
