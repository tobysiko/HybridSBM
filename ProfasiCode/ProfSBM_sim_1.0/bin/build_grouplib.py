#!/usr/bin/env python
import sys,re,string
libfile='residue_info.dat'
if len(sys.argv)>1: libfile=sys.argv[1]
libfile=open(libfile)
groupheader='GroupLib.hh'
groupsource='GroupLib.cc'
header=open(groupheader,'w')
source=open(groupsource,'w')
SPACES="                                                                      "
ind1=SPACES[:4]
ind2=SPACES[:8]
ind3=SPACES[:12]
ind4=SPACES[:16]
ind6=SPACES[:24]

print >> header, """\
#ifndef GroupLib_HH
#define GroupLib_HH
#include "GroupProps.hh"
#include "../Aux/profasi_io.hh"

//Automatically generated group library header file

namespace prf
{
"""
print >>source,"""\
#include "GroupLib.hh"

namespace prf
{

%snamespace Groups {
""" % (ind1)
lines=list()
OLC,chcode,startline,common_name,TLC,DESC=list(),list(),list(),list(),\
                                            list(),list()
linecount,groupcount=0,0
for line in libfile:
    if line.isspace(): continue
    lines.append(line)
    if re.match("^%",line):
        line= line.rstrip("\n")
        (marker,olc,charco,commn,tlc,desc)=re.split("\s+",line,5)
        map(lambda x,y:x.append(y),\
            (OLC,chcode,startline,common_name,TLC,DESC),\
            (olc,charco,linecount,commn,tlc,desc))
        groupcount+=1
    linecount+=1
libfile.close()
endline=startline[:]
endline.pop(0)
endline.append(len(lines))
print >>header,"""\
%(indnt)stypedef enum { %(olcs)s } OneLetterCode;

%(indnt)snamespace Groups 
%(indnt)s{
%(indnt2)sconst int max_olc=%(mxcodes)d;
%(indnt2)sextern GroupProps grp[max_olc];
%(indnt2)sextern char charcode[max_olc];
%(indnt2)sextern std::map<std::string,OneLetterCode> olcof;
%(indnt2)sbool checkGroup(std::string gg);
%(indnt2)svoid initGroups() ;
""" % {'indnt':ind1,'indnt2':ind2,'olcs':string.join(OLC,", "),'mxcodes':len(OLC)}
print >>source,"""\
%(indnt)sGroupProps grp[max_olc];
%(indnt)schar charcode[max_olc];
%(indnt)sstd::map<std::string,OneLetterCode> olcof;
%(indnt)sbool checkGroup(std::string gg)
%(indnt)s{
%(indnt2)sif (map2OLC(gg)==NONE) {
%(indnt3)sprf::cerr<<"Unknown group "<<gg<<"\\n";
%(indnt3)sreturn false;
%(indnt2)s}

%(indnt2)sreturn true;
%(indnt)s}

%(indnt)svoid initGroups()
%(indnt)s{
%(indnt2)sstatic bool already_initialized=false;

%(indnt2)sif (already_initialized) return;
"""  % {'indnt':ind2,'indnt2':ind3,'indnt3':ind4},

scdof=list()
labellist=list()
labellines=list()
atlst=False
for i in range(0,len(OLC)) :
    del scdof[:]
    del labellist[:]
    alllabels=""
    del labellines[:]
#    print "Line %d starts residue %s %s (%s) which is a %s" % (startline[i],OLC[i],common_name[i],TLC[i],DESC[i])
#    print "That residue ends in line %d" % (endline[i])
    for line in lines[startline[i]:endline[i]] :
        #print "For residue %d processing line %s" % (i,line),
        if re.match("^%",line) : continue
        elif re.match("atoms {",line) : atlst=True
        elif re.match("}",line) : atlst=False
        elif re.match("side_chain_dof",line) :
            (keywrd,dofid,at1,at2,at3,at4)=line.split()
            scdof.append(list((at1,at2,at3,at4)))
        elif atlst :
            line=line.rstrip()
            labellist.extend(re.split("\s+,\s+",line))
           # print line
    alllabels=string.join(labellist,"")
    alllabels=re.sub(r'"',r'',alllabels)
    alllabels=re.sub(r'\n',r'',alllabels)
    #print "Length of alllabels is ",len(alllabels)
    for k in range(0,len(alllabels)/40) :
        labellines.append(alllabels[40*k:40*(k+1)])
        #print "appended ",alllabels[40*k:40*(k+1)], " to label lines"
    labellines.append(alllabels[40*(len(alllabels)/40):])
    # print alllabels
    print >>source, '\n%s//%s' % (ind3,common_name[i])
    print >>source, '%sgrp[%s].CommonName("%s");\n' %(ind3,OLC[i],common_name[i])
    source.write(ind3+'grp['+OLC[i]+'].init(')
    nprev=''
    for k in range(0,len(labellines)-1) :
        print >>source,'%sstd::string("%s") +'% (nprev,labellines[k])
        nprev=SPACES[:len(OLC[i])+23]
    
    print >>source,'%sstd::string("%s"));\n' % (nprev,labellines[len(labellines)-1])
    #   grp[$curolc].init("$labelstr");
    if DESC[i].find("aminoacid")!=-1 :
        presentchar=chcode[i]
    else : 
        presentchar="X"
    print DESC[i],DESC[i].find("natural")
    form="""\
%(indnt)sgrp[%(curolc)s].set_type("%(curdes)s");

%(indnt)sgrp[%(curolc)s].TLC("%(curtlc)s");

%(indnt)scharcode[%(curolc)s]='%(curchar)s';

%(indnt)solcof["%(curtlc)s"]=%(curolc)s;

%(indnt)solcof["%(cn)s"]=%(curolc)s;
"""
    print >>source, form % {'indnt':ind3,'curolc':OLC[i], 'curdes':DESC[i], 'labelstr':alllabels, \
        'cn':common_name[i], 'curtlc':TLC[i],'curchar':presentchar},
    if DESC[i].find("natural")!=-1  and DESC[i].find("aminoacid")!=-1: 
        print >>source, '\n%solcof["%s"]=%s;' % (ind3,OLC[i],OLC[i]) 
    for k in range(0, len(scdof)) : 
        (at1,at2,at3,at4)=(scdof[k][0],scdof[k][1],scdof[k][2],scdof[k][3])
        print >>source,"\n%sgrp[%s].torsion_dof(%s,%s,%s,%s);" % (ind3,OLC[i],at1,at2,at3,at4)

print >>header,"""
%(indnt)sinline OneLetterCode mapTLC2OLC(std::string tlc) {return olcof[tlc];}

%(indnt)sinline OneLetterCode map2OLC(std::string nm) {return olcof[nm];}

%(indnt)sinline bool isAA(std::string st) {return grp[map2OLC(st)].isAA();}

%(indnt)sinline bool isEG(std::string st){return grp[map2OLC(st)].isEG();}

%(indnt)sinline std::string egname(std::string st)
%(indnt)s{
%(indnt)s    if (grp[map2OLC(st)].isEG()) return grp[map2OLC(st)].CommonName();
%(indnt)s    else return std::string("null");
%(indnt)s}

%(indnt)sinline Output & operator<<(Output & os, OneLetterCode olc) {return os<<((int) olc);}

%(indnt)sinline char mapOLC2Char(OneLetterCode olc) { return charcode[olc];}

%(indnt)sinline std::string mapOLC2Name(OneLetterCode olc) { return grp[olc].CommonName();}

%(indnt)sinline std::string mapOLC2TLC(OneLetterCode olc) { return grp[olc].TLC();}

%(indnt)sinline std::string SCALabel(OneLetterCode olc,int i) { return grp[olc].label(grp[olc].index(" CB ")+i);}

%(indnt)sinline OneLetterCode mapChar2OLC(char c) {return olcof[std::string(1,c)];}
""" % {'indnt':ind2}

print >> header, """\
%(indnt)s}
}

#endif
""" % {'indnt':ind1}
print >>source,"""\
%(i1)salready_initialized=true;
%(i2)s}
%(i3)s}
}
""" % {'i1':ind3,'i2':ind2,'i3':ind1}
header.close()
source.close()

