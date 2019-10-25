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

#include "prf_xml.hh"
#include "fileutils.hh"
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace prf_xml;
using namespace prf_utils;
typedef std::multimap<std::string, XML_Node*>::iterator mmiter;

XML_Node::XML_Node() : the_name("uninitialized"), mystate(OK), cmin(0), cmax(0) {}
XML_Node::~XML_Node()
{
    for (std::deque<XML_Node *>::iterator it=chldrn.begin();
         it!=chldrn.end(); ++it) {
        if (*it) delete(*it);
    }
}
XML_Node::XML_Node(const XML_Node &tg) : the_name(tg.the_name),
        the_value(tg.the_value), attr(tg.attr), childptr(tg.childptr),
        chldrn(tg.chldrn), mystate(tg.mystate), cmin(tg.cmin),
        cmax(tg.cmax) {}
XML_Node::XML_Node(std::string nm, std::string vl) : the_name(nm),
        the_value(vl), mystate(OK), cmin(0),cmax(0) {}

void XML_Node::make_clone_of(const XML_Node *tg)
{
    if (this!=tg) {
        the_name=tg->the_name;
        the_value=tg->the_value;
        mystate=tg->mystate;
        cmin=tg->cmin;
        cmax=tg->cmax;
        attr=tg->attr;
        for (size_t i=0;i<tg->n_children();++i) {
            XML_Node *newchild=new XML_Node();
            newchild->make_clone_of(tg->chldrn[i]);
            add_child_node(newchild);
        }
    }
}

XML_Node &XML_Node::operator=(const XML_Node &tg)
{
    if (this!=&tg) {
        the_name=tg.the_name;
        the_value=tg.the_value;
        attr=tg.attr;
        childptr=tg.childptr;
        chldrn=tg.chldrn;
        mystate=tg.mystate;
        cmin=tg.cmin;
        cmax=tg.cmax;
    }

    return *this;
}

XML_Node::XML_Node(std::list<prf_xml::Chunk>::iterator bg,
                   std::list<prf_xml::Chunk>::iterator nd)
{
    mystate=OK;
    std::string name1,name2;

    if (bg!=nd) {
        if (bg->type()==begin_tag) {
            name1=bg->tag_name();
            attr=bg->attribute_list();
        } else {
            prf::cerr<<"XML Node construction error: initial chunk \""
                     <<bg->show()<<"\" on line number "<<bg->line_number()
                     <<" is not a valid begin_tag.\n";
            mystate=ERROR;
        }

        if (nd->type()==end_tag) {
            name2=nd->tag_name();
        } else {
            prf::cerr<<"XML Node construction error: final chunk \""
                     <<nd->show()<<"\" on line number "<<bg->line_number()
                     <<" is not a valid end_tag.\n";
            mystate=ERROR;
        }

        if (mystate!=ERROR and name1!=name2) {
            prf::cerr<<"XML Node construction error: Begin tag "
                     <<bg->show()<<" on line number "<<bg->line_number()
                     <<" closed by mis-matched end tag "<<nd->show()
                     <<" on line number "<<nd->line_number()<<"\n";
            mystate=ERROR;
        }

        the_value="";

        for (std::list<prf_xml::Chunk>::iterator i=bg; i!=nd; ++i) {
            if (i->type()==text) the_value+=i->show();
        }

        cmin=bg->begin();
        cmax=nd->end();

        if (mystate!=ERROR) the_name=name1;
    } else {
        mystate=ERROR;
    }
}

void XML_Node::write(prf::Output &op, int indent_level)
{
    if (indent_level<0) indent_level=0;

    std::string ind(4*indent_level,' ');
    op<<ind<<the_name;

    if (n_attributes()!=0) {
        op<<"(";

        for (std::map<std::string,std::string>::const_iterator it=attr.begin();
             it!=attr.end(); ++it) {
            if (it!=attr.begin()) op<<", "; else op<<" ";

            op<<(it->first + "=\"" + it->second+"\"");
        }

        op<<")";
    }

    if (is_leaf_node()) {
        op<<" : \""<<the_value<<"\"\n";
        return;
    }

    op<<" {\n";
    ++indent_level;
    if (!the_value.empty()) {
        ind=std::string(4*indent_level,' ');
        op<<ind;

        for (size_t i=0; i<the_value.size(); ++i) {
            op<<the_value[i];

            if (the_value[i]=='\n') op<<ind;
        }

        op<<"\n";
    }
    if (n_children()>0) {
        for (std::deque<prf_xml::XML_Node *>::iterator it=chldrn.begin();
             it!=chldrn.end(); ++it) {
            (*it)->write(op,indent_level);
        }
    }

    --indent_level;
    ind=std::string(4*indent_level,' ');
    op<<ind<<"}\n";
}

std::string XML_Node::make_string() const
{
    std::string ans="<"+the_name;

    if (n_attributes()!=0) {
        for (std::map<std::string,std::string>::const_iterator it=attr.begin();
             it!=attr.end(); ++it) {
            ans+=(" "+ it->first + "=" + '"'+ it->second+'"');
        }
    }

    ans+=">\n";

    if (!the_value.empty()) ans+=the_value;

    for (std::deque<XML_Node *>::const_iterator it=chldrn.begin();
         it!=chldrn.end(); ++it) {
        if (*it) ans+=(*it)->make_string();
    }

    ans+="</"+the_name+">\n";
    return ans;
}

int XML_Node::add_child_node(XML_Node *ch)
{
    if (ch==NULL) return 0;

    if (child(ch->name())==ch) return 0;

    chldrn.push_back(ch);
    childptr.insert(std::pair<std::string,XML_Node *>(ch->name(),ch));
    return 1;
}

int XML_Node::remove_child_node(std::string nm)
{
    XML_Node *nd=child(nm);
    if (nd==NULL) return 0;
    mmiter mt;
    mt=childptr.find(nm);
    childptr.erase(mt);
    std::deque<XML_Node *>::iterator it=std::find(chldrn.begin(),
                                        chldrn.end(),nd);
    chldrn.erase(it);
    delete nd;
    return 1;
}

int XML_Node::disown_child_node(XML_Node *nd)
{
    if (nd==NULL) return 0;
    mmiter mt;
    std::pair<mmiter,mmiter> rnge;
    rnge=childptr.equal_range(nd->name());
    for (mt=rnge.first;mt!=rnge.second;++mt) {
        if (mt->second==nd) break;
    }
    if (mt!=rnge.second) childptr.erase(mt);
    std::deque<XML_Node*>::iterator it;
    it=std::find(chldrn.begin(),chldrn.end(),nd);
    if (it!=chldrn.end()) chldrn.erase(it); else return 0;
    return 1;
}

int XML_Node::remove_child_node(XML_Node *nd)
{
    if (disown_child_node(nd)) {
        delete nd;
        return 1;
    } else return 0;
}

void XML_Node::clear_child_nodes()
{
    childptr.clear();
    for (size_t i=0;i<chldrn.size();++i) {
        delete chldrn[i];
        chldrn[i]=NULL;
    }
    chldrn.clear();
}

int XML_Node::add_child_node(std::string nm, std::string vl)
{
    XML_Node *chld=new XML_Node(nm,vl);
    return add_child_node(chld);
}

int XML_Node::interpret_formatted_data()
{
    std::deque<XML_Node *> newchldrn;
    for (size_t i=0;i<chldrn.size();) {
        XML_Node *nd=chldrn[i];
        if (nd==NULL or nd->name()!="formatted_data") {
            if (nd!=NULL) nd->interpret_formatted_data();
            ++i;
            continue;
        }
        XML_Node *fmt=nd->child("format");
        std::string tmplstr=fmt->make_string();
        XML_Node *dt=nd->child("data");
        XML_Node *imp=nd->child("import_data");
        if (fmt==NULL or (dt==NULL && imp==NULL)) {
            ++i;
            continue;
        }
        std::string line,fileimport="",nodeimport="";
        if (dt!=NULL) nodeimport=dt->value();
        if (imp!=NULL) {
            prf_utils::get_file_contents(imp->value(),fileimport);
        }
        std::istringstream ist;
        nodeimport+=fileimport;
        ist.str(nodeimport);
        while (getline(ist,line)) {
            line=prf_utils::trim_str(line);
            if (!line.empty()) {
                std::deque<std::string> parts;
                prf_utils::split(line,parts);
                std::string newtmpl=tmplstr;
                std::istringstream inp;
                std::ostringstream out;
                inp.str(newtmpl);
                char c;
                while (inp.get(c)) {
                    if (c!='$') out<<c;
                    else {
                        size_t fldno=0;
                        inp>>fldno;
                        if (fldno==0) out<<line;
                        else if ((fldno-1)<parts.size()) out<<parts[fldno-1];
                    }
                }
                XML_Mini xml;
                std::string processeddata=out.str();
                xml.parse_data(processeddata);
                XML_Node *rcd=xml.build_tree();
                rcd->set_name(fmt->attribute("name"));
                rcd->attributes().erase("name");
                newchldrn.push_back(rcd);
            }
        }
        if (remove_child_node(nd)==0) ++i;
    }

    for (size_t i=0;i<newchldrn.size();++i) add_child_node(newchldrn[i]);
    newchldrn.clear();

    return 1;
}

Chunk::Chunk() : i1(0),i2(0), line_no(1), mytype(text) {}
Chunk::~Chunk() {}
Chunk::Chunk(const Chunk &c) : the_text(c.the_text), i1(c.i1), i2(c.i2),
        line_no(c.line_no), mytype(c.mytype) {}
Chunk::Chunk(std::string txt, size_t bg, size_t nd) : the_text(txt),
        i1(bg), i2(nd), line_no(1), mytype(text) {}
Chunk &Chunk::operator=(const Chunk &c)
{
    if (this!=&c) {
        the_text=c.the_text;
        i1=c.i1;
        i2=c.i2;
        line_no=c.line_no;
        mytype=c.mytype;
    }

    return *this;
}
std::string Chunk::info()
{
    std::string ans("Chunk: ");

    switch (mytype) {
        case bad: ans+="ILL FORMED "; break;
        case text: ans+="Text "; break;
        case begin_tag: ans+="Begin tag "; break;
        case end_tag: ans+="End tag "; break;
        case special_tag: ans+="Special tag "; break;
        default: break;
    };

    ans+="Limits ";

    char lim[50];

    sprintf(lim, "(%u, %u)", (unsigned int) i1, (unsigned int) i2);

    ans+=std::string(lim);

    ans+=("Content:"+the_text);

    return ans;
}

std::string Chunk::tag_name()
{
    return tagname;
}

std::map<std::string, std::string> Chunk::attribute_list()
{
    if (mytype!=begin_tag) {
        prf::cerr<<"Failed to get attributes. "<<the_text
                 <<" is not a valid begin tag.\n";
    }

    return atl;
}

chunk_type Chunk::determine_type()
{
    the_text=trim_str(the_text);
    size_t sz=the_text.size();
    mytype=bad;

    if (the_text[0]=='>') return mytype;

    if (the_text[0]=='<' && the_text[sz-1]!='>') return mytype;

    if (the_text[0]=='<') {
        if (the_text[1]=='?' || the_text[1]=='!') mytype=special_tag;
        else if (the_text[1]=='/') {
            std::string tagtext=the_text.substr(2,the_text.size()-3);

            //Check if it is an end tag
            if (isalpha(the_text[2]) and valid_name(tagtext)) {
                mytype=end_tag;
                tagname=tagtext;
            }
        } else {
            //Check if it is a begin tag
            std::string tagtext,atributes;
            size_t i=1;

            while (i<the_text.size() and
                   (isalpha(the_text[i]) or isdigit(the_text[i]) or
                    the_text[i]=='_')) tagtext+=the_text[i++];

            std::string attributes=trim_str(the_text.substr(i,the_text.size()-i-1));

            if ((!tagtext.empty()) and isalpha(tagtext[0]) and
                syntax_check_attributes(attributes, atl)) {
                mytype=begin_tag;
                tagname=tagtext;
            }
        }
    } else if (syntax_check_text(the_text)) mytype=text;
    else prf::cerr<<"Invalid chunk \""<<the_text<<"\" at line number "
                      <<line_no<<"\n";

    return mytype;
}

bool prf_xml::syntax_check_text(std::string txt)
{
    bool good_text=true;

    for (size_t i=0; good_text && i<txt.size(); ++i) {
        //Strictly speaking '>' is allowed in XML. But we don't allow it
        //for simplicity.
        if (txt[i]=='<' or txt[i]=='>') good_text=false;

        if (txt[i]=='&') {
            std::string er; //entity reference

            for (size_t j=i+1; j<txt.size(); ++j) {
                er+=txt[j];

                if (not isalpha(txt[j])) break;
            }

            if (er.empty() or
                (er!="lt;" && er!="gt;" && er!="amp;"
                 && er!="apos;" && er!="quot;")) good_text=false;
        }
    }

    return good_text;
}

bool prf_xml::valid_name(std::string nm)
{
    //1. Must begin with a normal character
    //2. Can contain only alphanumeric characters and '_'
    nm=trim_str(nm);
//     std::cout<<"Checking if \""<<nm<<"\" is a valid name\n";
    bool alnum=(not nm.empty()) and isalpha(nm[0]);

    for (size_t i=1; alnum&&i<nm.size(); ++i) {
        alnum=isalpha(nm[i]) or isdigit(nm[i]) or nm[i]=='_';
    }

    return alnum;
}

bool prf_xml::get_atr_name(std::string att, size_t &ipos, std::string &nm)
{
    nm="";
    while (ipos<att.size() and att[ipos]!='=') nm+=att[ipos++];
    nm=trim_str(nm);
    return (ipos<att.size() and att[ipos]=='=');
}

bool prf_xml::get_atr_text(std::string att, size_t &ipos, std::string &tx)
{
    tx="";
    while (ipos<att.size() and isspace(att[ipos])) ++ipos;
    if (not(ipos<att.size() and (att[ipos]=='\"' or att[ipos]=='\'')))
        return false;
    char openqot=att[ipos];
    bool closed=false;
    while ((++ipos)<att.size() and not closed) {
        if (att[ipos]==openqot) closed=true;
        else tx+=att[ipos];
    }
    if (not closed) tx.clear();
    return closed;
}

bool prf_xml::syntax_check_attributes(std::string att,
                                      std::map<std::string, std::string> &atl)
{
    bool syntax_ok=true;
    att=trim_str(att);
    atl.clear();

//     std::cout<<"Checking if \""<<att<<"\" forms a good attributes spec.\n";
    if (att.empty()) return syntax_ok;

    std::string atname,attext;
    size_t ipos=0;
    do {
        bool foundeq=get_atr_name(att,ipos,atname);
        if (atname.empty()) {
            syntax_ok=(not foundeq);
            break;
        }
        if (not (syntax_ok=valid_name(atname))) break;
        ++ipos;
        syntax_ok=get_atr_text(att,ipos,attext)
                    and syntax_check_text(attext);
        if (syntax_ok) atl[atname]=attext;
    } while (syntax_ok);

    if (!syntax_ok) atl.clear();
    return syntax_ok;
}

XML_Mini::XML_Mini() {}
XML_Mini::~XML_Mini() {}
void XML_Mini::clear()
{
    data.clear();
}

size_t XML_Mini::read_file(std::string fl)
{
    if (!TestFile_r(fl)) return 0;
    std::string filcont;
    prf_utils::get_file_contents(fl,filcont);
    return parse_data(filcont);
}

size_t XML_Mini::parse_data(std::string &dat) {
    enum {text_input_mode,tag_input_mode} input_mode;
    char c;
    std::string tmp="";
    size_t bg=0, nd=0,lno=1,cloc=0;
    data.clear();
    input_mode=text_input_mode;

    for (cloc=0;cloc<dat.size();++cloc) {
        c=dat[cloc];
        if (input_mode==text_input_mode) {
            if (c=='<') {
                if (!trim_str(tmp).empty()) data.push_back(Chunk(tmp,bg,nd));

                tmp.clear();
                bg=nd;
                input_mode=tag_input_mode;
            }

            tmp+=c;
            nd++;
        } else {
            tmp+=c;
            nd++;

            if (c=='>') {
                Chunk chu(tmp,bg,nd);
                chu.line_number(lno);
                data.push_back(chu);
                tmp.clear();
                bg=nd;
                input_mode=text_input_mode;
            }
        }

        if (c=='\n') ++lno;
    }

    tmp=trim_str(tmp);

    if (!tmp.empty()) data.push_back(Chunk(tmp,bg,nd));

    for (std::list<prf_xml::Chunk>::iterator it=data.begin();
         it!=data.end(); ++it) {
        it->determine_type();
    }

    return data.size();
}

XML_Node *XML_Mini::build_tree()
{
    XML_Node *the_root=NULL;
    std::list<Chunk> special_tags, faulty_tags;
    std::list<Chunk>::iterator it,jt;
    it=jt=data.begin();

    while (jt!=data.end()) {
        ++jt;

        if (it->type()==special_tag) special_tags.splice(special_tags.end(),data,it);
        else if (it->type()==bad) faulty_tags.splice(faulty_tags.end(),data,it);

        it=jt;
    }

    if (data.empty()) return the_root;

    if (!faulty_tags.empty()) {
        prf::cerr<<"XML Error! "<<faulty_tags.size()<<" faulty tags!\n";
        int i=0;

        for (it=faulty_tags.begin(); it!=faulty_tags.end(); ++it,++i) {
            prf::cerr<<i<<": "<<it->show()<<"\n";
        }

        return the_root;
    }

    size_t stagnation_meter=0;
    std::list<XML_Node *> nodes;
    state_type parser_state=OK;

    while (!data.empty()) {
        size_t prev_size=data.size();
        std::list<prf_xml::Chunk>::iterator last_begin_tag, t,u;
        t=u=data.begin();
        last_begin_tag=data.end();

        while (u!=data.end()) {
            ++u;

            if (t->type()==begin_tag) last_begin_tag=t;
            else if (t->type()==end_tag && last_begin_tag!=data.end()) {
                XML_Node *nnd= new XML_Node(last_begin_tag,t);

                if (nnd->state()==ERROR) {
                    delete nnd;
                    parser_state=ERROR;
                } else {
                    std::list<XML_Node *>::iterator ni,nj;
                    ni=nj=nodes.begin();

                    while (nj!=nodes.end()) {
                        ++nj;

                        if (nnd->encloses(**ni)) {
                            nnd->add_child_node(*ni);
                            nodes.erase(ni);
                        }

                        ni=nj;
                    }

                    nodes.push_back(nnd);
                }

                data.erase(last_begin_tag,u);
                last_begin_tag=data.end();
            }

            t=u;
        }

        if (data.size()!=prev_size) {
            stagnation_meter=0;
            prev_size=data.size();
        } else ++stagnation_meter;

        if (parser_state==ERROR or stagnation_meter>5) break;
    }

    if (nodes.size()!=1) {
        prf::cerr<<"XML Error! Found "<<nodes.size()
                 <<" nodes at the top level. There must be one and "
                 <<" only one root element. \n";

        if (nodes.size()>1) {
            prf::cerr<<"This could happen if the file actually contains "
                     <<"more than one top level elements, or if the parsing "
                     <<"is aborted due to errors before the whole file is "
                     <<"processed.\n";
        }
    } else {
        the_root=nodes.front();
    }

    return the_root;
}

XML_Node *XML_Node::child(std::string nm)
{
    std::multimap<std::string, XML_Node *>::iterator it=childptr.begin();
    XML_Node *ans=NULL;

    while (it!=childptr.end() && it->first!=nm) ++it;

    if (it!=childptr.end()) ans=it->second;

    return ans;
}

XML_Node *prf_xml::get_xml_tree(std::string xml_file_name)
{
    XML_Mini xml;
    xml.read_file(xml_file_name);
    return xml.build_tree();
}

XML_Node *prf_xml::make_xml_tree(std::string xml_node_text)
{
    XML_Mini xml;
    xml.parse_data(xml_node_text);
    return xml.build_tree();
}
