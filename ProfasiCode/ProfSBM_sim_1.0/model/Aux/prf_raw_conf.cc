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

#include "prf_raw_conf.hh"
#include "fileutils.hh"
#include <deque>
#include <sstream>

using namespace prf_traj;
using namespace prf_utils;
using std::ios;

prf_raw_conf::prf_raw_conf()
{
    byte_order_of_file="Unknown";
    datafl="uninitiailzed_conf_file";
    infofl="uninitiailzed_conf_file.info";
    mymap=NULL;
    buf=NULL;
    signature.clear();
    sizes_host.resize(7,0);
    sizes_host[CHAR]=sizeof(char);
    sizes_host[INT]=sizeof(int);
    sizes_host[LONG]=sizeof(long);
    sizes_host[FLOAT]=sizeof(float);
    sizes_host[DOUBLE]=sizeof(double);
    sizes_host[LONG_DOUBLE]=sizeof(long double);
    sizes_host[UNKNOWN]=0;
    sizes_file=sizes_host;
    init_xml_map();
    awake=false;
}

prf_raw_conf::prf_raw_conf(const prf_raw_conf &cf)
{
    datafl=cf.datafl;
    infofl=cf.infofl;
    headertxt=cf.headertxt;
    byte_order_of_file=cf.byte_order_of_file;
    blocksz=cf.blocksz;
    nblk=cf.nblk;
    ncrd=cf.ncrd;
    header_offset=cf.header_offset;
    nbytes=cf.nbytes;
    sizes_file=cf.sizes_file;
    sizes_host=cf.sizes_host;
    if (blocksz!=0) buf=new char[blocksz]; else buf=NULL;
    signature=cf.signature;
    if (cf.mymap!=NULL) {
        mymap=new prf_xml::XML_Node();
        mymap->make_clone_of(cf.mymap);
    }
    awake=initialized=false;
}

prf_raw_conf &prf_raw_conf::operator=(const prf_raw_conf &cf)
{
    if (this!=&cf) {
        datafl=cf.datafl;
        infofl=cf.infofl;
        headertxt=cf.headertxt;
        byte_order_of_file=cf.byte_order_of_file;
        blocksz=cf.blocksz;
        nblk=cf.nblk;
        ncrd=cf.ncrd;
        header_offset=cf.header_offset;
        nbytes=cf.nbytes;
        sizes_file=cf.sizes_file;
        sizes_host=cf.sizes_host;
        if (blocksz!=0) buf=new char[blocksz]; else buf=NULL;
        signature=cf.signature;
        if (cf.mymap!=NULL) {
            if (mymap!=NULL) delete mymap;
            mymap=new prf_xml::XML_Node();
            mymap->make_clone_of(cf.mymap);
        }
        awake=initialized=false;
    }
    return *this;
}
prf_raw_conf::~prf_raw_conf()
{
    if (myfl.is_open()) myfl.close();
    if (buf) delete buf;
    if (mymap) delete mymap;
}

void prf_raw_conf::set_file(std::string st)
{
    if (st!=datafl or (!initialized)) {
        initialized=false;
        unsigned loc=st.find_last_of('/');
        if (loc>=st.size()) infofl="conf.info";
        else infofl=st.substr(0,loc)+"/conf.info";
    }
    datafl=st;
}

bool prf_raw_conf::init()
{
    if (initialized) return wake_up();
    init_xml_map();
//    prf::cout<<"Opening file "<<datafl.c_str()<<"\n";
    myfl.open(datafl.c_str());
    if ((!myfl.good()) or (!get_header()) or (!parse_header())) {
        if (myfl.is_open()) myfl.close();
        return awake=false;
    }
    if (blocksz==0) {
        prf::cerr<<"prf_raw_conf> Found block size of 0 although header "
                <<"has apparently been processed!\n";
        return false;
    }
    if (nbytes<(header_offset+blocksz)) {
        prf::cerr<<"Not enough data for even one complete record in "
                <<datafl<<"\n";
        return false;
    }
    nblk=(nbytes-header_offset)/blocksz;
    if ((nbytes-header_offset)%blocksz!=0) {
        prf::cerr<<"The file "<<datafl<<" contains corrupt data.\n";
        prf::cerr<<"Block size is "<<blocksz<<" where as number of bytes "
                <<"in the binary data is "
                <<(nbytes-header_offset)<<"\n";
    }
    buf=new char[blocksz];
    awake=true;
    return initialized=true;
}

bool prf_raw_conf::sleep()
{
    if (awake) {
        if (myfl.is_open()) myfl.close();
        awake=false;
    }
    return true;
}

bool prf_raw_conf::wake_up()
{
    if (!awake) {
        if (!initialized) {
            return init();
        } else {
            if (!myfl.is_open()) myfl.open(datafl.c_str());
            if (!myfl.good()) {
                prf::cerr<<"prf_raw_conf> Failing to wake up as the file "
                        <<datafl<<" could not be opened. \n";
                return awake=false;
            } else return awake=true;
        }
    }
    return true;
}

void prf_raw_conf::init_xml_map()
{
    if (mymap) delete mymap;
    mymap=new prf_xml::XML_Node();
    mymap->set_name("profasi_conf_map");
    mymap->add_child_node("population","");
}

bool prf_raw_conf::get_block(unsigned iblk, TrajSnapshot &snp)
{
    if (iblk>=nblk) return false;
    myfl.seekg(header_offset+iblk*blocksz);
    if (myfl.fail()) return false;
    myfl.read(buf,blocksz);
    if (myfl.fail()) return false;
    bin.data(buf,blocksz);
    snp.set_n_coordinates(ncrd);
    snp.fill(bin);
    return true;
}

bool prf_raw_conf::get_header()
{
    return read_inline_header() or read_info_file();
}

bool prf_raw_conf::read_inline_header()
{
    myfl.seekg(0,ios::end);
    nbytes=myfl.tellg();
    myfl.seekg(0,ios::beg);

    bool header_exists=true;
    if (nbytes<21) header_exists=false;
    else {
        char headertag[20];
        myfl.get(headertag,20);
        if (std::string(headertag)!="PROFASI_CONF_HEADER") {
            header_exists=false;
            myfl.seekg(0,ios::beg);
            headertxt="";
        } else {
            headertxt=std::string(headertag)+"\n";
        }
    }

    if (header_exists) {
        std::string line;
        while (getline(myfl,line)) {
            if (!line.empty()) headertxt+=(line+"\n");
            if (prf_utils::rm_space(line)=="END_PROFASI_CONF_HEADER") break;
        }
    }
    header_offset=myfl.tellg();
    return header_exists;
}

bool prf_raw_conf::read_info_file()
{
    return prf_utils::get_file_contents(infofl,headertxt)!=0;
}

bool prf_raw_conf::parse_header()
{
    std::deque<std::string> lines;
    std::istringstream ssin(headertxt);
    std::string line;
    while (getline(ssin,line)) lines.push_back(line);
    int curch=0;
    unsigned iPs=0,iPe=0;
    for (size_t i=0;i<lines.size();++i) {
        line=lines[i];
        line=prf_utils::rm_space(line);
        if (line.empty()) continue;
        std::deque<std::string> parts;
        prf_utils::split<std::deque<std::string> >(line,parts);
        std::string keyword=parts.front();
        if (keyword=="conf_length") {
            blocksz=atoi(parts[1].c_str());
            continue;
        } else if (keyword=="byte_order") {
            byte_order_of_file=parts[1];
        } else if (keyword=="sequence") {
            prf_xml::XML_Node *newch=new prf_xml::XML_Node();
            newch->set_name("protein");
            std::ostringstream sout;
            sout<<curch++;
            newch->set_attribute("id",sout.str());
            newch->add_child_node("sequence",line.substr(8));
            mymap->child("population")->add_child_node(newch);
        } else if (keyword.size()>10 and keyword.substr(0,7)=="peptide" and
                   keyword.substr(keyword.size()-5)=="start") {
            if (curch==0) iPs=atoi(parts[1].c_str());
        } else if (keyword.size()>10 and keyword.substr(0,7)=="peptide" and
                   keyword.substr(keyword.size()-3)=="end") {
            iPe=atoi(parts[1].c_str());
        } else if (keyword=="box_length") {
            mymap->add_child_node("box_length",parts[1]);
        } else if (keyword=="ntmp") {
            mymap->add_child_node("num_temperatures",parts[1]);
        } else if (keyword=="temperature_array") {
            std::ostringstream ost;
            for (size_t it=1;it<parts.size();++it) {
                ost<<parts[it];
                if (it!=(parts.size()-1)) ost<<" ";
            }
            mymap->add_child_node("temperature_array",ost.str());
        } else if (keyword=="sizeof_char")
            sizes_file[CHAR]=atoi(parts[1].c_str());
        else if (keyword=="sizeof_int")
            sizes_file[INT]=atoi(parts[1].c_str());
        else if (keyword=="sizeof_long")
            sizes_file[LONG]=atoi(parts[1].c_str());
        else if (keyword=="sizeof_float")
            sizes_file[FLOAT]=atoi(parts[1].c_str());
        else if (keyword=="sizeof_double")
            sizes_file[DOUBLE]=atoi(parts[1].c_str());
        else if (keyword=="sizeof_long_double")
            sizes_file[LONG_DOUBLE]=atoi(parts[1].c_str());
        int pos2=1;
        VarType tp=UNKNOWN;
        if (keyword=="unsigned") keyword+=(" "+parts[pos2++]);
        if ((tp=prf_utils::make_type(keyword))!=UNKNOWN) {
            int noftyp=atoi(parts[pos2].c_str());
            for (int j=0;j<noftyp;++j) signature.push_back(tp);
        }
    }
    bin.set_data_attribs(byte_order_of_file,sizes_file);
    ncrd=(iPe-iPs)/sizes_file[DOUBLE];
    std::ostringstream sout;
    sout<<curch;
    mymap->child("population")->add_child_node("num_chains",sout.str());
    unsigned long checksum=0;
    for (size_t i=0;i<signature.size();++i) checksum+=sizes_file[signature[i]];
    if (checksum!=blocksz) {
        prf::cerr<<"Sum of individual sizes = "<<checksum<<"\n";
        prf::cerr<<"Block size given in file = "<<blocksz<<"\n";
        return false;
    }
    check_sizes();
    return true;
}

bool prf_raw_conf::check_sizes()
{
    for (size_t i=0;i<sizes_file.size();++i) {
        if (sizes_file[i]>sizes_host[i]) {
            if (std::find(signature.begin(),signature.end(),(VarType) i)
                !=signature.end()) {
                prf::cerr<<"This machine stores "<<make_type_str(i)
                        <<" in "<<sizes_host[i]<<" bytes, where as "
                        <<"the machine where the conf file was written stores"
                        <<" them in "<<sizes_file[i]<<" bytes. If the values "
                        <<"read from the conf files are too big for this "
                        <<"machine, they will be truncated. This could lead "
                        <<"to problems.\n";
            }
        }
    }
    return true;
}


