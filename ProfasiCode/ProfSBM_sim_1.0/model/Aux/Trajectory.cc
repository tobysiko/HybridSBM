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

#include "Trajectory.hh"
#include <algorithm>
#include <fstream>

using namespace prf_traj;

TrajIterator::TrajIterator() : idx(0), isgood(false),
needsrefresh(true), mytraj(0) {}

TrajIterator::TrajIterator(const TrajIterator &t) : idx(t.idx), mytraj(t.mytraj)
{
    isgood=t.isgood;
    needsrefresh=t.needsrefresh;
    snpsh=t.snpsh;
}

TrajIterator::TrajIterator(Trajectory *tr, unsigned long wh) : idx(wh),
mytraj(tr)
{
    isgood=false;
    needsrefresh=true;
}

TrajIterator::~TrajIterator() {}

TrajIterator & TrajIterator::operator=(const TrajIterator &t)
                                      {
    if (this!=&t) {
        idx=t.idx;
        mytraj=t.mytraj;
        isgood=t.isgood;
        needsrefresh=t.needsrefresh;
        snpsh=t.snpsh;
    }
    return *this;
}

bool TrajIterator::operator==(const TrajIterator &t)
{
    return idx==t.idx and mytraj==t.mytraj;
}

bool TrajIterator::operator !=(const TrajIterator &t)
{
    return idx!=t.idx or mytraj!=t.mytraj;
}

bool TrajIterator::operator<(const TrajIterator &t)
{
    return idx<t.idx and mytraj==t.mytraj;
}

bool TrajIterator::operator>(const TrajIterator &t)
{
    return idx>t.idx and mytraj==t.mytraj;
}

bool TrajIterator::operator<=(const TrajIterator &t)
{
    return idx<=t.idx and mytraj==t.mytraj;
}

bool TrajIterator::operator>=(const TrajIterator &t)
{
    return idx>=t.idx and mytraj==t.mytraj;
}

int TrajIterator::operator+=(int n)
                            {
    isgood=false;
    needsrefresh=true;
    idx+=n;
    if (idx<mytraj->end().idx) return 0;
    unsigned long bkp=idx;
    idx=mytraj->end().idx;
    return bkp-idx+1;
}

TrajIterator & TrajIterator::operator++()
{
    ++idx;
    isgood=false;
    needsrefresh=true;
    return *this;
}

TrajIterator TrajIterator::operator++(int)
{
    TrajIterator tmp(*this);
    ++idx;
    isgood=false;
    needsrefresh=true;
    return tmp;
}

TrajIterator &TrajIterator::operator--()
{
    --idx;
    isgood=false;
    needsrefresh=true;
    return *this;
}

TrajIterator TrajIterator::operator--(int)
{
    TrajIterator tmp(*this);
    --idx;
    isgood=false;
    needsrefresh=true;
    return *this;
}

TrajSnapshot & TrajIterator::operator*()
{
    refresh();
    return snpsh;
}

TrajSnapshot * TrajIterator::operator->()
{
    refresh();
    return &snpsh;
}

void TrajIterator::refresh()
{
    if (needsrefresh) {
        isgood=mytraj->focus(idx,snpsh);
        if (isgood) needsrefresh=false;
    }
}

Trajectory::Trajectory() : curblk(0), curseg(0), cycmin(0), cycmax(0),
bgptr(this,0), ndptr(this,0), myfile("traj")
{
    mymap=new prf_xml::XML_Node("MCTrajectory","");
}

Trajectory::~Trajectory()
{
    for (size_t i=0;i<segment.size();++i) {
        if (segment[i]) {
            delete segment[i];
            segment[i]=NULL;
        }
    }
    if (mymap) delete mymap;
}

void Trajectory::save(std::string flnm)
{
    save_up_to(max_cycle(),flnm);
}

void Trajectory::save_up_to(unsigned long icyc, std::string flnm)
{
    prf::Output trj(flnm.c_str(),"w");

    std::list<std::string> opath,ipath;
    prf_utils::make_rel_path(flnm,opath);

    std::string obasename=flnm;
    if (!opath.empty()) {
        obasename=opath.back();
        opath.pop_back();
    }

    trj<<"#_PROFASI_TRAJECTORY_DESCRIPTION_FILE\n";
    trj<<"# LAYOUT: SEGMENT_FILE START_CYCLE INTERVAL N_BLOCKS\n";
    for (size_t i=0;i<segment.size();++i) {
        TrajSeg *seg=segment[i];
        if (seg->min_cycle()<=icyc) {
            prf_utils::make_rel_path(seg->get_filename(),ipath);
            std::list<std::string>::iterator ot,it;
            ot=opath.begin();
            it=ipath.begin();
            while (ot!=opath.end() and it!=ipath.end() and *ot==*it) {
                ++ot;
                ++it;
            }
            std::string relname;
            for (std::list<std::string>::iterator kt=ot;kt!=opath.end();++kt) {
                relname=relname+"../";
            }
            for (std::list<std::string>::iterator kt=it;kt!=ipath.end();++kt) {
                if (kt!=it) relname=relname+"/";
                relname=relname+(*kt);
            }
            trj<<relname<<" "<<seg->min_cycle()<<" "
                <<seg->interval()<<" "<<seg->n_blocks()<<"\n";
        }
    }
    trj.close();
}

int Trajectory::parse(std::string trajfile)
{
    if (prf_utils::TestFile_r(trajfile)==0) return 0;
    std::ifstream myfl(trajfile.c_str());
    myfl.seekg(0,std::ios::end);
    unsigned long nbytes=myfl.tellg();
    myfl.seekg(0,std::ios::beg);

    char headertag[38];
    if (nbytes>=38) myfl.get(headertag,38);
    if (nbytes<38 or
        std::string(headertag)!="#_PROFASI_TRAJECTORY_DESCRIPTION_FILE") {
        prf::cerr<<"Trajectory> File "<<trajfile<<" lacks a valid header "
                <<"for a PROFASI trajectory description file.\n";
        return -1;
    }
    std::string refdir,line;
    size_t loc=trajfile.find_last_of('/');
    if (loc<trajfile.size()) refdir=trajfile.substr(0,loc);
    else refdir=".";
    std::deque<std::string> lines;
    while (getline(myfl,line)) lines.push_back(line);
    myfl.close();
    int nseg=append_list(refdir, lines);
    if (nseg>0) return 1; else return 0;
}

int Trajectory::append_list(std::string refdir, std::deque<std::string> &lines)
{
    int nseg=0;
    for (size_t i=0;i<lines.size();++i) {
        if (lines[i].empty() or lines[i][0]=='#') continue;
        std::deque<std::string> parts;
        prf_utils::split(lines[i],parts);
        if (parts.size()==1) {
            nseg+=append_segment(refdir+"/"+parts[0]);
        } else if (!parts.empty()) {
            unsigned long strt=strtoul(parts[1].c_str(),NULL,10);
            unsigned intv=abs(atoi(parts[2].c_str()));
            unsigned nblks=abs(atoi(parts[3].c_str()));
            nseg+=append_segment(refdir+"/"+parts[0],strt,intv,nblks);
        }
    }
    return nseg;
}

int Trajectory::init()
{
    if (segment.empty()) return 0;
    return auto_adjust();
}

int Trajectory::auto_adjust()
{
    if (segment.empty()) return 0;
    std::vector<std::pair<unsigned long,unsigned> > segorder(segment.size());
    for (unsigned i=0;i<segorder.size();++i) {
        segorder[i].first=segment[i]->min_cycle();
        segorder[i].second=i;
    }
    std::sort(segorder.begin(),segorder.end());
    std::list<TrajSeg *> tmpseg;
    for (unsigned i=0;i<segorder.size();++i) {
        tmpseg.push_back(segment[segorder[i].second]);
    }
    std::list<TrajSeg *>::iterator i1,i2;
    i2=i1=tmpseg.begin();
    if (i2!=tmpseg.end()) ++i2;
    while (i2!=tmpseg.end()) {
        unsigned long cyclim=(*i2)->min_cycle();
        unsigned long mymincyc=(*i1)->min_cycle();
        unsigned long mymaxcyc=(*i1)->max_cycle();

        if (cyclim<=mymincyc) {
            delete (*i1);
            tmpseg.erase(i1);
        } else if (mymaxcyc>=cyclim) {
            (*i1)->truncate_in_t_range(mymincyc,cyclim-1);
        }
        i1=i2++;
    }

    if (tmpseg.empty()) return 0;
    segment.resize(tmpseg.size(),NULL);
    std::copy(tmpseg.begin(),tmpseg.end(),segment.begin());

    unsigned iblk=0;
    for (size_t i=0;i<segment.size();++i) {
        segment[i]->set_block_id_offset(iblk);
        iblk+=segment[i]->n_blocks();
    }
    ndptr=iterator(this,iblk);
    cycmin=segment.front()->min_cycle();
    cycmax=segment.back()->max_cycle();
    if (!segment.empty()) {
        prf_xml::XML_Node *segnode=segment[0]->xml_map();
        if (segnode!=NULL) {
            for (size_t i=0;i<segnode->n_children();++i) {
                std::string nm=segnode->child(i)->name();
                if (nm=="population" or nm=="box_length" or
                    nm=="ntmp" or nm=="temperature_array") {
                    mymap->remove_child_node(nm);
                    prf_xml::XML_Node *newnode=new prf_xml::XML_Node();
                    newnode->make_clone_of(segnode->child(i));
                    mymap->add_child_node(newnode);
                }
            }
        }
    }
    return 1;
}

bool Trajectory::focus(unsigned iblk, prf_traj::TrajSnapshot &snpsh)
{
    unsigned newseg=seg_of_block(iblk);
    if (newseg>=segment.size()) return false;
    if (newseg!=curseg) {
        segment[curseg]->sleep();
        if (!segment[newseg]->wake_up()) return false;
        curseg=newseg;
    }
    if (segment[curseg]->get_block(iblk,snpsh)) {
        curblk=iblk;
        return true;
    }
    return false;
}

unsigned Trajectory::seg_of_block(unsigned blkid)
{
    if (blkid>=segment[curseg]->first_block() &&
        blkid<=segment[curseg]->last_block()) return curseg;
    if (blkid>segment[curseg]->last_block()) {
        unsigned newseg=curseg+1;
        while (newseg<segment.size() && segment[newseg]->last_block()<blkid)
            ++newseg;
        return newseg;
    } else {
        unsigned newseg=curseg;
        while (newseg>0) {
            --newseg;
            if (segment[newseg]->first_block()<=blkid) break;
        }
        if (segment[newseg]->first_block()<=blkid) return newseg;
        else return segment.size();
    }
}

Trajectory::iterator Trajectory::find(unsigned long mct)
{
    unsigned newseg=seg_of_time(mct),newblk=0;
    if (newseg<segment.size()) {
        if (segment[newseg]->calc_block_id(mct,newblk))
            return iterator(this, newblk);
    }
    return ndptr;
}

unsigned Trajectory::seg_of_time(unsigned long mct)
{
    if (mct>=segment[curseg]->min_cycle()) {
        unsigned newseg=curseg+1;
        while (newseg<segment.size() && segment[newseg]->min_cycle()<=mct)
            ++newseg;
        return newseg-1;
    } else {
        unsigned newseg=curseg;
        while (newseg>0) {
            --newseg;
            if (segment[newseg]->min_cycle()<=mct) break;
        }
        if (segment[newseg]->min_cycle()<=mct) return newseg;
        else return segment.size();
    }
}

int Trajectory::append_segment(std::string sgfl)
{
    TrajSeg *seg=new TrajSeg();
    seg->set_file(sgfl);
    if (!seg->init()) {
        delete seg;
        return 0;
    }
    if (curseg<segment.size()) {
        segment[curseg]->sleep();
    }
    segment.push_back(seg);
    curseg=segment.size()-1;
    return 1;
}

int Trajectory::append_segment(std::string sgfl, unsigned long st,
                               unsigned long iv, unsigned long nbl)
{
    TrajSeg *seg;
    if (segment.empty()) {
        seg=new TrajSeg();
        seg->set_file(sgfl);
        if (seg->init()) {
            curseg=0;
            segment.push_back(seg);
        } else {
            delete seg;
            return 0;
        }
    } else {
        seg=new TrajSeg(*segment[0]);
        seg->set_file(sgfl);
        seg->set_block_props(st,iv,nbl);
        segment.push_back(seg);
    }
    return 1;
}
