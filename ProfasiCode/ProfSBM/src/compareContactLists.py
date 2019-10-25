import glob,sys,os.path,math,ast
from matplotlib import pylab as p
from matplotlib import markers as ma

AA_3_to_1 = {
					'ALA':'A',
					'ARG':'R',
					'ASN':'N',
					'ASP':'D',
					'CYS':'C',
					'GLN':'Q',
					'GLU':'E',
					'GLY':'G',
					'HIS':'H',
					'ILE':'I',
					'LEU':'L',
					'LYS':'K',
					'MET':'M',
					'PHE':'F',
					'PRO':'P',
					'SER':'S',
					'THR':'T',
					'TRP':'W',
					'TYR':'Y',
					'VAL':'V'}
def mergeStrings(s1,s2):
	s=""
	for i in xrange(min(len(s1),len(s2))):
		if s1[i]==s2[i]:
			s += s1[i]
		else:
			s += "x"
	return s
def writeContactHisto(n, histo, filename):
	f = open(filename.replace("/","_"),"w")
	f.write(str(n) + "\n")
	f.write(str(histo) + "\n")
	f.close()
	print "wrote histogram to",filename
	return True
def readContactHisto(hname):
	f = open(hname)
	n = int(f.readline())
	histo = ast.literal_eval(f.read())
	print "read histogram from",hname
	return histo,n

def getResidueDict(atoms):
	rd = {}
	
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in rd:
			rd[rnum].append(atoms[a])
		else:
			rd[rnum]=[atoms[a]]
	
	
	return rd
def atomDist(a1,a2,atoms):
	d = math.sqrt( math.pow(atoms[a1]["x"] - atoms[a2]["x"], 2)
					+ math.pow(atoms[a1]["y"] - atoms[a2]["y"], 2)
					+ math.pow(atoms[a1]["z"] - atoms[a2]["z"], 2)
					)
	#print d
	return d
def parseATOMline(line, returnDict=True): 
	anum = int(line[6:11])
	aname = line[12:16].strip()
	altloc = line[16]
	rname = line[17:20].strip()
	chain = line[21]
	rnum = int(line[22:26])
	insert = line[26]
	x = float(line[30:38])
	y = float(line[38:46])
	z = float(line[46:54])
	occupancy = float(line[54:60])
	bfactor = float(line[60:66])
	element = line[76:78].strip()
	charge = line[78:80].strip()
	
	if returnDict:
		return {"anum":anum, "aname":aname, "altloc":altloc, "rname":rname, "chain":chain, "rnum":rnum, "insert":insert, "x":x, "y":y, "z":z, "occupancy":occupancy, "bfactor":bfactor, "element":element, "charge":charge}
	else:
		return (anum, aname, altloc, rname, chain, rnum, insert, x, y, z, occupancy, bfactor, element, charge)

def atomsFromPDB(filename):
	atoms = {}
	for line in open(filename).readlines():
		if line[:4]=="ATOM":
			a = parseATOMline(line)
			atoms[a["anum"]] = a
	return atoms
def getResidueContacts(atoms, cutoff, chainDist, verbose=False, res2SS = None, includeNatoms=False):
	#from Bio.PDB import *
	
	#structure = PDBParser.get_structure(PDBParser(),pdbfilename,pdbfilename)

	cutoff = cutoff

	chainsep = chainDist

	#model = structure[0]

	#reslist = Selection.unfold_entities(structure, 'R')

	contacts = []
	residues = getResidueDict(atoms)
	n=len(residues)
	#print n,"residues"
	
	for i in xrange(n):
		ri = residues[i+1]
		for j in xrange(n):
			rj = residues[j+1]
		
			if i<j-chainsep:
				mindist = 10000
				mdpair = None
				cadist = None
				natoms = 0
				for ai in ri:
					if not (ai["aname"][0]=='H' or (ai["aname"][0] in ['1','2','3','4'] and ai["aname"][1]=='H') or ai["aname"] in ["CA","C","O","N"]): 
						#if ai["rname"]=="GLY": print ai["aname"],ai["rname"],ai["rnum"],"accepted"
						for aj in rj:
							if not (aj["aname"][0]=='H' or (aj["aname"][0] in ['1','2','3','4'] and aj["aname"][1]=='H') or aj["aname"] in ["CA","C","O","N"]): 
								
								d = atomDist(ai["anum"],aj["anum"],atoms)
					
								if cadist == None:
									if ai["aname"] == 'CA' and aj["aname"]== 'CA':
										cadist = d
								if d < mindist:
									mindist=d
									mdpair = (ai["anum"],aj["anum"])
									natoms += 1
							else:
								pass
					else:
						#if ai["rname"]=="GLY": print ai["aname"],ai["rname"],ai["rnum"]
						pass
				if mindist <= cutoff:
					#a1 = getAtomIDforResidueContact(i+1,'CA',atoms)
					#a2 = getAtomIDforResidueContact(j+1,'CA',atoms)
					
					if includeNatoms:
						contacts.append( (i+1,j+1,natoms,mindist) )
					else:
						contacts.append( (i+1,j+1) )
					#print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	
	
	return contacts

def contactHisto(pdbfiles, cutoff, chainDist): 
	
	contactHistoAll = {}
	
	for pdbfile in pdbfiles:
		atoms = atomsFromPDB(pdbfile)
		contacts = getResidueContacts(atoms, cutoff, chainDist, verbose=False, res2SS = None)
		for c in contacts:
			
			if c in contactHistoAll:
				contactHistoAll[c] += 1
			else:
				contactHistoAll[c] = 1
	
	for h in contactHistoAll:
		contactHistoAll[h] /= float(len(pdbfiles))
	
	return contactHistoAll

def getResidueDistance(structure, r1, r2, c1, c2):
	mindist = 10e10
	minpair = None
	
	atoms1 = structure[0][c1][r1].get_list()
	atoms2 = structure[0][c2][r2].get_list()
	
	
	for i in xrange(len(atoms1)):
		a1 = atoms1[i]
		for j in xrange(len(atoms2)):
			a2 = atoms2[j]
			
			if i<j:
				d = a1-a2
				if d < mindist:
					mindist = d
					minpair = (a1,a2)
	return mindist,minpair

def readXMLContacts(filename):
	contacts = {}
	start = False
	for line in open(filename).readlines():
		if "<data>" in line:
			start = True
			continue
		elif "</data>" in line:
			break
		
		if start:
			tmp = line.strip().split()
			if len(tmp) < 2: continue
			tmp1 = tmp[0].split("/")
			tmp2 = tmp[1].split("/")
			
			r1 = int(tmp1[1])+1
			r2 = int(tmp2[1])+1
			contacts[(r1,r2)] = tmp
	return contacts

def getPDBfileList(name,filt=None,readlinks=False):
	flist=[]
	if os.path.isdir(name):
		if filt==None:
			flist = glob.glob(name+"/*.pdb")
		else:
			if readlinks:
				flist = glob.glob(name+"/*.pdb")
				tmp = []
				for f in flist:
					l = os.path.basename(os.readlink(f))
					if filt in l:
						tmp.append(f)
				flist = tmp
				#print flist
			else:
				flist = glob.glob(name+"/%s*.pdb"%filt)
	elif os.path.isfile(name) and ".pdb" in name:
		flist = [name]
	else:
		print "ERROR: Unknown pdb file or list of pdb files: ",name
		sys.exit(1)
	
	
	
	return flist
def getLength(atoms):
	maxrnum = 0
	for a in atoms:
		if atoms[a]["rnum"]>maxrnum:
			maxrnum = atoms[a]["rnum"]
	return maxrnum
def getSequence(atoms):
	rdic = getResidueDict(atoms)
	s=[]
	for r in sorted(rdic.keys()):
		s.append(AA_3_to_1[rdic[r][0]['rname'].upper()])
			
	return s
###################################################################



inpdb1 = sys.argv[1]


if len(sys.argv) == 3:
	inpdb2 = sys.argv[2]
	compareByName = False
else:
	inpdb2 = inpdb1
	compareByName = True

filterlabels = ["GA","GB"]

alwaysOverwrite = False

readlinks = True

distcut = 6.0
chainsep = 3

chain1 = "A"
chain2 = "A"

nativeContactFile1 = "superGA.pdb_FMULTIGAUSS_RESIDUE_c6.0_n3_w0.5_R0.0_CA_cons1.0_xFalse_95contacts_E65.0.xml"
nativeContactFile2 = "superGB.pdb_FMULTIGAUSS_RESIDUE_c6.0_n3_w0.5_R0.0_CA_cons1.0_xFalse_137contacts_E68.0.xml"

native1 = readXMLContacts(nativeContactFile1)
native2 = readXMLContacts(nativeContactFile2)

print "num. native contacts:",len(native1),len(native2)

if compareByName:
	pdbfiles1 = getPDBfileList(inpdb1,filt=filterlabels[0],readlinks=readlinks)
	pdbfiles2 = getPDBfileList(inpdb2,filt=filterlabels[1],readlinks=readlinks)
	inpdb1 = filterlabels[0]+inpdb1
	inpdb2 = filterlabels[1]+inpdb2
else:
	pdbfiles1 = getPDBfileList(inpdb1,readlinks=readlinks)
	pdbfiles2 = getPDBfileList(inpdb2,readlinks=readlinks)

print

print "num. pdb files:",len(pdbfiles1),len(pdbfiles2)

atoms1 = atomsFromPDB(pdbfiles1[0])
atoms2 = atomsFromPDB(pdbfiles2[0])

sequence1 = getSequence(atoms1)
sequence2 = getSequence(atoms2)
print "".join(sequence1)
print "".join(sequence2)

hname1 = inpdb1 + "_contacthisto.dat"
if not os.path.exists(hname1) or  alwaysOverwrite:
	histo1 = contactHisto(pdbfiles1, distcut, chainsep)
	
	
	nres1 = getLength(atoms1)
	
	assert writeContactHisto(nres1, histo1, hname1)
else:
	
	histo1,nres1 = readContactHisto(hname1)



hname2 = inpdb2 + "_contacthisto.dat"
if not os.path.exists(hname2) or alwaysOverwrite:
	histo2 = contactHisto(pdbfiles2, distcut, chainsep)
	
	nres2 = getLength(atoms2)
	assert writeContactHisto(nres2, histo2, hname2)
	
else:
	histo2,nres2 = readContactHisto(hname2)

print "num. residues:",nres1,nres2
assert nres1==nres2


allcontactpairs = {}

for h in histo1.keys():
	if not h in allcontactpairs:
		allcontactpairs[h]=True
for h in histo2.keys():
	if not h in allcontactpairs:
		allcontactpairs[h]=True


print "num. contacts in histogram:",len(histo1),len(histo2)


datafile = open(mergeStrings(inpdb1,inpdb2).replace("/","_")+"_contactDiff.csv",'w')


#print "pos1,pos2,co,p1,p2,dp,nat1?,nat2?"
datafile.write("pos1,pos2,co,p1,p2,dp,nat1?,nat2?\n")

rows = []

p.figure(dpi=100)
for c in allcontactpairs:
	if c in histo1:
		p1 = histo1[c]
	else:
		p1 = 0.0
	if c in histo2:
		p2 = histo2[c]
	else:
		p2 = 0.0
	
	res1 = c[0]
	res2 = c[1]
	co = abs(c[0] - c[1])
	dp= p2-p1
	maxp = max(p1,p2)
	nat1 = int(c in native1)
	nat2 = int(c in native2)
	
	row = [res1, res2, co, p1, p2, dp, nat1, nat2]
	rows.append(row)
	#print ",".join([str(i) for i in row])
	
	
	#mark = ma.MarkerStyle()
	#mark.set_marker("o")
	#mark.set_fillstyle('none')
	
	#print mark.get_fillstyle(), mark.get_marker()
	if nat1 and not nat2:
		col = "blue"
	elif nat2 and not nat1:
		col = "red"
	elif nat1 and nat2:
		col = "magenta"
	else:
		col = "black"
	p.scatter([dp], [co], s=[maxp*500.0], color=col,edgecolor=col,facecolor=(0,0,0,0))
	
	import random
	if abs(dp) >= 0.02 and maxp >= 0.2:
		#m=math.log(co)-math.log(10)
		m=co-10
		u= m+20
		l=m-40
		
		
		if dp > 0:
			
			offset = (random.randint(20,30),random.randint(l,u))
			seq=sequence2
		elif dp < 0:
			
			offset = (random.randint(-70,-50),random.randint(l,u))
			seq=sequence1
		else:
			xoffset = 0.0
			seq=sequence2
		p.annotate("(%s%i,%s%i)"%(seq[res1-1],res1,seq[res2-1],res2),(dp,co),color=col,size=8,textcoords='offset points',xytext=offset,arrowprops={'width':0.1,'frac':0.1,'headwidth':1,'shrink':0.05})
		#p.text(dp,co,"(%s%i,%s%i)"%(seq[res1-1],res1,seq[res2-1],res2),color=col)

rows.sort(key=lambda x: x[5], reverse=True)

for r in rows:
	datafile.write(",".join([str(i) for i in r])+"\n")

datafile.close()
#ymin, ymax = p.ylim()
p.annotate("A",(10,10),xytext=(0.1,0.9),textcoords='axes fraction',color='blue')
p.annotate("B",(10,10),xytext=(0.9,0.9),textcoords='axes fraction',color='red')
p.ylim(chainsep-1, nres1+10)
p.grid()
#p.legend(loc=0)
p.xlabel("P(B) - P(A)")
p.ylabel("contact order")
p.title("A: %s\nB: %s"%(inpdb1,inpdb2))
p.yscale('log')
p.yticks([4,5,6,7,8,9,10,20,30,40,50,60],[4,5,6,7,8,9,10,20,30,40,50,60])

outfname = mergeStrings(inpdb1,inpdb2)
p.savefig(outfname.replace("/","_")+"_contact_diff.png",dpi=200)
#p.show()



