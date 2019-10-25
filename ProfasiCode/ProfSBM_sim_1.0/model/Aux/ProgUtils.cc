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

#include "ProgUtils.hh"
#include "profasi_io.hh"
#include "fileutils.hh"

using namespace prf;

using std::list;
using std::vector;
using std::deque;
using std::string;
using std::map;

namespace prf_utils
{
    OptBase::OptBase()
    {
        myname="unknown";
        longname="unknown";
        shhlp="";
        hlp="";
        given=false;
    }

    OptBase::~OptBase() {}

    OptBase::OptBase(std::string longnm, std::string shortnm, std::string shorthlp)
    {
        longname=longnm;
        myname=shortnm;
        shhlp=shorthlp;
        given=false;
    }

    OptBase::OptBase(const OptBase &ob)
    {
        longname=ob.longname;
        myname=ob.myname;
        shhlp=ob.shhlp;
        hlp=ob.hlp;
        given=ob.given;
    }

    OptBase & OptBase::operator=(const OptBase &ob)
    {
        if (this!=&ob) {
            longname=ob.longname;
            myname=ob.myname;
            shhlp=ob.shhlp;
            hlp=ob.hlp;
            given=ob.given;
        }
        return *this;
    }

    ProgSwitch::ProgSwitch() : OptBase(), val(true) {}

    ProgSwitch::~ProgSwitch() {}

    ProgSwitch::ProgSwitch(std::string longnm, std::string shortnm,
                           bool defstate, std::string shorthlp) : OptBase(longnm,shortnm,shorthlp)
    {
        val=defstate;
    }

    ProgSwitch::ProgSwitch(const ProgSwitch &ps) : OptBase(ps)
    {
        val=ps.val;
    }

    ProgSwitch & ProgSwitch::operator=(const ProgSwitch &ps)
    {
        if (this!=&ps) {
            OptBase::operator=(ps);
            val=ps.val;
        }
        return *this;
    }

    ProgOpt::ProgOpt() : OptBase()
    {
        na=0;
        vals.clear();
    }

    ProgOpt::ProgOpt(string lng, string srt, int i, string sth) : OptBase(lng,srt,sth)
    {
        na=i;
        vals.resize(i,"");
    }

    ProgOpt::ProgOpt(const ProgOpt &po) : OptBase(po)
    {
        na=po.na;
        vals=po.vals;
    }

    ProgOpt & ProgOpt::operator=(const ProgOpt &po)
    {
        if (this!=&po) {
            OptBase::operator =(po);
            na=po.na;
            vals=po.vals;
        }
        return *this;
    }

    ProgOpt::~ProgOpt() {}

    ProgArgs::ProgArgs()
    {
        settingsfile="settings.cnf";
        rank=0;
        opts.clear();
        remn.clear();
        available.clear();
        option("settings_file","st",1,
               "(Name of settings file, if it is not \"settings.cnf\")");
        settings_use_type(0);
    }

    ProgArgs::~ProgArgs() {}

    void ProgArgs::set_availability(string opnm, bool avl)
    {
        deque<ProgOpt>::iterator it=opts.begin();

        for (;it!=opts.end();++it)
            if (it->Name()==opnm ||it->LongName()==opnm) break;

        if (it!=opts.end()) {
            available[it->Name()]=avl;
            available[it->LongName()]=avl;
        }

        deque<ProgSwitch>::iterator st=swtchs.begin();

        for (;st!=swtchs.end();++st)
            if (st->Name()==opnm ||st->LongName()==opnm) break;

        if (st!=swtchs.end()) {
            available[st->Name()]=avl;
            available[st->LongName()]=avl;
        }
    }

    void ProgArgs::disable(std::string opnm)
    {
        set_availability(opnm,false);
    }

    void ProgArgs::enable(string opnm)
    {
        set_availability(opnm,true);
    }

    void ProgArgs::option(string longname,string shortname, int nargus,
                          string helptext)
    {
        if (nargus>0) opts.push_back(ProgOpt(longname,shortname,nargus, helptext));
        else swtchs.push_back(ProgSwitch(longname,shortname,true,helptext));
        available[shortname]=true;
        available[longname]=true;
    }

    void ProgArgs::new_switch(std::string longnm,std::string shortnm,
                              bool defstate, std::string shorthlp)
    {
        swtchs.push_back(ProgSwitch(longnm,shortnm,defstate,shorthlp));
        available[shortnm]=true;
        available[longnm]=true;
    }

    string ProgArgs::option_arr(string optname, int vali)
    {
        deque<ProgOpt>::iterator it;

        for (it=opts.begin();it!=opts.end();++it) {
            if (it->Name()==optname or it->LongName()==optname) break;
        }

        if (it!=opts.end()) {
            if (vali>=it->max_args()) vali=it->max_args()-1;

            return it->value(vali);
        }

        return string("");
    }

    string ProgArgs::option(string optname)
    {
        return option_arr(optname,0);
    }

    bool ProgArgs::option_given(string optname)
    {
        deque<ProgOpt>::iterator it;

        for (it=opts.begin();it!=opts.end();++it) {
            if (it->Name()==optname or it->LongName()==optname) break;
        }

        if (it!=opts.end()) {
            return it->state();
        } else {
            return false;
        }
    }

    bool ProgArgs::switch_given(std::string swnm)
    {
        std::deque<ProgSwitch>::iterator it;

        for (it=swtchs.begin();it!=swtchs.end();++it) {
            if (it->Name()==swnm or it->LongName()==swnm) break;
        }

        if (it!=swtchs.end()) {
            return it->state();
        } else return false;
    }

    bool ProgArgs::state_of_switch(std::string swnm)
    {
        std::deque<ProgSwitch>::iterator it;

        for (it=swtchs.begin();it!=swtchs.end();++it) {
            if (it->Name()==swnm or it->LongName()==swnm) break;
        }

        if (it!=swtchs.end()) {
            return it->on();
        } else return false;
    }

    string ProgArgs::spare_args(int i)
    {
        if (i<(int)remn.size()) return remn[i];
        else prf::cerr<<"Unused parameter number "<<i<<" asked for. \n"
            <<"But only "<<remn.size()<<" parameters are unused after "
            <<"parsing command line options.\n";

        return string("");
    }

    void ProgArgs::write_available()
    {
        if (opts.size()>0) prf::cout<<"Options taking one or more arguments...\n";
        for (deque<ProgOpt>::iterator it=opts.begin();it!=opts.end();++it) {
            if (available[it->Name()])
                prf::cout<<"\t-"<<it->Name()<<"  or --"<<it->LongName()<<"\t: "
                <<it->max_args()<<" argument"<<(it->max_args()==1?". ":"s. ")
                <<it->short_help()<<"\n";
        }
        if (swtchs.size()>0) prf::cout<<"Program switches, which can be given as on or off...\n";
        for (deque<ProgSwitch>::iterator it=swtchs.begin();it!=swtchs.end();++it) {
            if (available[it->Name()])
                prf::cout<<"\t-"<<it->Name()<<"  or --"<<it->LongName()<<"\t: "
                <<(it->on()?"true":"false")<<it->short_help()<<"\n";
        }
    }

    void ProgArgs::analyze(int argc,char *argv[])
    {
        list<string> cmdln;

        for (int i=1;i<argc;++i) {
            cmdln.push_back(string(argv[i]));
        }

        for (list<string>::iterator i=cmdln.begin();i!=cmdln.end();) {
            string wrd=*i;
            string chopped;
            bool uselongname=false,optfound=false,flgfound=false,flgnegated=false;

            if ((!is_number(wrd)) && wrd[0]=='-') {
                int chopbegin=1;

                if (wrd[1]=='-') {
                    uselongname=true; chopbegin=2;
                }

                chopped=string(wrd,chopbegin);

                deque<ProgOpt>::iterator it;
                deque<ProgSwitch>::iterator st;

                if (uselongname) {
                    for (it=opts.begin();it!=opts.end();++it) {
                        if (it->LongName()==chopped) {
                            optfound=true;
                            break;
                        }
                    }
                    if (!optfound) {
                        for (st=swtchs.begin();
                        st!=swtchs.end(); ++st) {
                            if (st->LongName()==chopped) {
                                flgfound=true;
                                break;
                            } else if (("no-"+st->LongName())==chopped) {
                                flgfound=flgnegated=true;
                                break;
                            }
                        }
                    }
                } else {
                    for (it=opts.begin();it!=opts.end();++it) {
                        if (it->Name()==chopped) {
                            optfound=true;
                            break;
                        }
                    }
                    if (!optfound) {
                        for (st=swtchs.begin();
                        st!=swtchs.end(); ++st) {
                            if (st->Name()==chopped) {
                                flgfound=true;
                                break;
                            } else if (("no-"+st->Name())==chopped) {
                                flgfound=flgnegated=true;
                                break;
                            }
                        }
                    }
                }

                list<string>::iterator j=i;++j;

                if (optfound) {
                    int nwargs=0;

                    while (j!=cmdln.end() &&
                           (is_number(*j) || (*j)[0]!='-') &&
                           nwargs<it->max_args()) {
                        it->value(nwargs++,*j++);
                    }

                    if (nwargs<it->max_args()) {
                        prf::cerr<<"Insufficient number of arguments passed to"
                        <<" option "<<it->Name()<<". \n"
                        <<"Option ignored.\n";
                    } else if (available[it->Name()]) {
                        it->set_given();
                        std::string tmpstr=it->LongName()+" ";
                        if (it->max_args()==0) tmpstr+="on";
                        else
                            for (int j=0;j<it->max_args();++j) tmpstr+=(it->value(j)+" ");
                        cmdlni.push_back(InstructionString(tmpstr));
                    }
                } else if (flgfound) {
                    if (available[st->Name()]) {
                        st->set_given();
                        if (flgnegated) {
                            st->turn_off();
                            cmdlni.push_back(InstructionString(st->LongName()+" off"));
                        } else {
                            st->turn_on();
                            cmdlni.push_back(InstructionString(st->LongName()+" on"));
                        }
                    }
                }

                i=cmdln.erase(i,j);
            } else ++i;
        }

        int ii=0;

        remn.resize(cmdln.size());

        for (list<string>::iterator i=cmdln.begin();i!=cmdln.end();++i) {
            remn[ii++]=*i;
        }
    }

    int ProgArgs::get_primary_settings()
    {
#ifdef MPI
        int bufsz=0;
        std::string sbuf;
        if (rank==0) bufsz=get_file_contents(settingsfile,sbuf);
        MPI::Comm::Bcast(&bufsz,1,MPI::INT,0);
        if (bufsz==0) return 0;
        char *buf=new char[bufsz+1];
        for (size_t i=0;i<sbuf.size();++i) buf[i]=sbuf[i];
        buf[bufsz]='\0';
        MPI::Comm::Bcast(&buf,bufsz+1,MPI::CHAR,0);
        istringstream ssin(buf);
        get_commands(ssin,cmds);
#else
        if (STestFile_r(settingsfile.c_str())==0) return 0;
        get_commands(settingsfile,cmds);
#endif
        return 1;
    }

    int ProgArgs::get_secondary_settings()
    {
        bool secsettings=false;
        for (std::list<InstructionString>::iterator it=cmds.begin();
        it!=cmds.end();++it) {
            if (it->head()=="secondary_settings" and it->tail().str()=="on") {
                secsettings=true;
                break;
            }
        }
        char sf[50];
        if (!secsettings) return 0;
        sprintf(sf,"n%d/",rank);
        std::string secsf=sf+settingsfile;
        if (TestFile_r(secsf)==0) return 0;
        std::list<InstructionString> tmpstngs;
        get_commands(secsf,tmpstngs);
        cmds.splice(cmds.end(),tmpstngs);
        return 1;
    }

    int ProgArgs::get_settings()
    {
        get_primary_settings();
        get_secondary_settings();
        for (std::list<InstructionString>::iterator it=cmds.begin(); it!=cmds.end();++it) {
            if (it->head()=="for_rank") {
                int trnk=atoi(it->part(1).c_str());
                if (rank==trnk) {
                    (*it)=(it->tail().tail());
                }
            }
        }
        return 1;
    }

    void ProgArgs::init_options(int argc, char *argv[])
    {
        analyze(argc,argv);
        if (option_given("settings_file")) settingsfile=option("settings_file");
        if (settings_use==1 or
            (settings_use==2 and option_given("st"))) get_settings();
        add_cmd_line_instructions();
    }

    void ProgArgs::add_cmd_line_instructions()
    {
        cmds.splice(cmds.end(),cmdlni);
    }

    bool ProgArgs::option_available(std::string optname)
    {
        bool ans=false;
        for (size_t i=0;i<opts.size();++i) {
            if (((opts[i].Name()==optname) or
                 (opts[i].LongName()==optname)) and opts[i].state()) {
                ans=true;
                break;
            }
        }
        return ans;
    }
}


