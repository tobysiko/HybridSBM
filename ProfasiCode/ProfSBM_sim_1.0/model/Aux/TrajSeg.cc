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

#include "TrajSeg.hh"

using namespace prf_traj;

TrajSeg::TrajSeg() : mncyc(0), intv(1000), mxcyc(0), fmncyc(0), fmxcyc(0),
blkidoff(0), nblkfile(0), nblk(0),frstusedblk(0), awake(false)
{
    mymap=NULL;
    data=NULL;
}

TrajSeg::~TrajSeg()
{
    if (data) delete data;
    if (mymap) delete mymap;
}

TrajSeg::TrajSeg(const TrajSeg &t)
{
    mncyc=t.mncyc;
    mxcyc=t.mxcyc;
    intv=t.intv;
    fmncyc=t.fmncyc;
    fmxcyc=t.fmxcyc;
    blkidoff=t.blkidoff;
    nblkfile=t.nblkfile;
    nblk=t.nblk;
    frstusedblk=t.frstusedblk;
    if (t.mymap!=NULL) {
        mymap=new prf_xml::XML_Node();
        mymap->make_clone_of(t.mymap);
    } else mymap=NULL;
    if (t.data!=NULL) data=new prf_raw_conf(*t.data);
    else data=NULL;
    awake=false;
}

bool TrajSeg::get_block(unsigned iblk, TrajSnapshot &snp)
{
    if (iblk<blkidoff or iblk>=(blkidoff+nblk)) return false;
    return get_file_block(iblk-blkidoff+frstusedblk, snp);
}

bool TrajSeg::get_block_at_mct(unsigned long mct, TrajSnapshot &snp)
{
    unsigned tmpblk=0;
    if (!calc_block_id(mct,tmpblk)) return false;
    return get_block(tmpblk,snp);
}

bool TrajSeg::calc_block_id(unsigned long mct, unsigned &blkid)
{
    if (mct<mncyc) return false;
    if (mct>(mxcyc+intv-1)) return false;
    unsigned tmpblk=(mct-mncyc)/intv;
    if (tmpblk<nblk) {
        blkid=blkidoff+tmpblk;
        return true;
    }
    return false;
}

void TrajSeg::truncate_in_block_range(unsigned bl1, unsigned bl2)
{
    if (bl1>bl2) std::swap(bl1,bl2);
    if (bl1>nblkfile) bl1=nblkfile;
    frstusedblk=bl1;
    if (bl2>nblkfile) bl2=nblkfile;
    nblk=bl2-bl1+1;
    mncyc=fmncyc+frstusedblk*intv;
    mxcyc=mncyc;
    if (nblk>0) {
        mxcyc+=(nblk-1)*intv;
    }
}

void TrajSeg::truncate_in_t_range(unsigned long t0, unsigned long t1)
{
    if (t0>t1) std::swap(t0,t1);
    if (t0<fmncyc) frstusedblk=0; else {
        frstusedblk=(t0-fmncyc)/intv;
    }
    unsigned tmplblk=nblkfile-1;
    if (t1<=fmxcyc) {
        tmplblk=(t1-fmncyc)/intv;
    }
    truncate_in_block_range(frstusedblk,tmplblk);
}

std::string TrajSeg::determine_conf_type(std::string st)
{
    std::string ans="empty";
    if (!st.empty()) ans="PROFASI_RAW_CONF";
    return ans;
}

void TrajSeg::set_file(std::string st)
{
    std::string filetype=determine_conf_type(st);
    if (filetype=="PROFASI_RAW_CONF") {
        if (data) delete data;
        data=new prf_raw_conf();
    }
    data->set_file(st);
}

std::string TrajSeg::get_filename()
{
    if (data) return data->get_filename();
    else return "null_segment";
}

void TrajSeg::set_block_props(unsigned long cyc0, unsigned long intvl,
                              unsigned blks)
{
    fmncyc=mncyc=cyc0;
    nblk=nblkfile=blks;
    intv=intvl;
    fmxcyc=mxcyc=mncyc+intv*(nblk-1);
}

int TrajSeg::init()
{
    if (!data->init()) return 0;
    if (mymap) delete mymap;
    mymap=new prf_xml::XML_Node();
    mymap->make_clone_of(data->xml_map());
    nblk=nblkfile=data->n_blocks();
    if (nblk==0) return 0;
    blkidoff=frstusedblk=0;
    prf_traj::TrajSnapshot s1,s2;
    data->get_block(0,s1);
    data->get_block(nblk-1,s2);
    mncyc=fmncyc=s1.MC_time();
    mxcyc=fmxcyc=s2.MC_time();
    if (nblk==1) intv=10000000;
    else intv=(fmxcyc-fmncyc)/(nblk-1);
    return 1;
}

bool TrajSeg::get_file_block(unsigned iblk, TrajSnapshot &snp)
{
    if (iblk>=nblkfile) return false;
    if (not awake) wake_up();
    return data->get_block(iblk,snp);
}

bool TrajSeg::sleep()
{
    return data->sleep();
}

bool TrajSeg::wake_up()
{
    return data->wake_up();
}
