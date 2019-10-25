import sys,math

def getAtomIDforResidueContact(resid,atype,atoms):
	
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype:
			
			return a
	
	
	return None

def getResidueContacts(atoms, cutoff, chainDist, verbose=False, res2SS = None):
	#from Bio.PDB import *
	
	#structure = PDBParser.get_structure(PDBParser(),pdbfilename,pdbfilename)

	cutoff = cutoff

	chainsep = chainDist

	#model = structure[0]

	#reslist = Selection.unfold_entities(structure, 'R')

	contacts = []
	residues = getResidueDict(atoms)
	offset = min(residues.keys())-1
	n=len(residues)
	print n,"residues"
	
	for i in xrange(n):
		ri = residues[i+1+offset]
		for j in xrange(n):
			rj = residues[j+1+offset]
		
			if i<j-chainsep:
				mindist = 10000
				mdpair = None
				cadist = None
				for ai in ri:
					if 'H' in ai["aname"]: 
						#print ai.get_id()
						continue
					for aj in rj:
						if 'H' in aj["aname"]: 
							continue
						d = atomDist(ai["anum"],aj["anum"],atoms)
					
						if cadist == None:
							if ai["aname"] == 'CA' and aj["aname"]== 'CA':
								cadist = d
						if d < mindist:
							mindist=d
							mdpair = (ai["anum"],aj["anum"])
				if mindist <= cutoff:
					a1 = getAtomIDforResidueContact(i+1,'CA',atoms)
					a2 = getAtomIDforResidueContact(j+1,'CA',atoms)
					
					contacts.append( [i+1+offset,j+1+offset,cadist] )
					#print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	
	
	return contacts

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

def plotCompContactMap(n, contacts1, contacts2, title = "", fname="contactmap.png"):
	
	import matplotlib.pyplot as p
	import matplotlib.cm as cm
	
	cdict1 = {}
	cdict2 = {}
	
	c2inc1 = 0
	for c in contacts1:
		r1 = int(c[0])
		r2 = int(c[1])
		#d = float(c[2])
		
		cdict1[(r1,r2)] = 1.0
		
	for c in contacts2:
		r1 = int(c[0])
		r2 = int(c[1])
		#d = float(c[2])
		
		if (r1,r2) in cdict1 or (r2,r1) in cdict1:
			c2inc1 += 1
			cdict2[(r1,r2)] = 1.0
		else:
			cdict2[(r1,r2)] = 0.2
	
	print c2inc1, "contacts of pdb1 found in pdb2"
	matrix = [[0.0 for j in xrange(n+1)] for i in xrange(n+1) ]
	
	for i in xrange(n+1):
		
		for j in xrange(n+1):
			
			
			
			if (i,j) in cdict1:
				matrix[i][j] = cdict1[(i,j)]
			
			
			if (i,j) in cdict2:
				matrix[j][i] = cdict2[(i,j)]
				
		
	
	title = title+"\n%i shared contacts"%c2inc1#"%s; cutoff=%s; csep=%i; n=%i\n%s"%(inPDBFilename,str(round(cutoff,2)),chainSep,len(contacts),contactType)

	imax = p.matshow(matrix,cmap=cm.binary,aspect='equal', origin='lower')
	p.plot([1,n],[1,n],color='k')
	p.title(title)
	#cbar = p.colorbar()
	p.xlabel("residue position PDB1")
	p.ylabel("residue position PDB2")
	p.xlim(0.5,n+0.5)
	p.ylim(0.5,n+0.5)
	p.grid(color="0.5", linestyle=':', linewidth=1.0)
	#cbar.set_label("CA dist")
	#cbar.set_label("CA dist. range")
	#imax.get_axes().set_axisbelow(True)
	p.savefig(fname,dpi=600)
	p.show()

def atomsFromPDB(filename):
	atoms = {}
	for line in open(filename).readlines():
		if line[:4]=="ATOM":
			a = parseATOMline(line)
			atoms[a["anum"]] = a
	return atoms

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

def getLength(atoms):
	maxrnum = 0
	for a in atoms:
		if atoms[a]["rnum"]>maxrnum:
			maxrnum = atoms[a]["rnum"]
	return maxrnum

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
			contacts[(r1,r2)] = tmp[3].split(",")[0]
			
	return contacts

pdb1 = sys.argv[1]
pdb2 = sys.argv[2]



chainsep = 3
distcut = 6

if pdb1[-4:]==".pdb":
	atoms1 = atomsFromPDB(pdb1)
	
	contacts1 = getResidueContacts(atoms1, distcut, chainsep)
elif pdb1[-4:]==".xml":
	contacts1 = readXMLContacts(pdb1)

atoms2 = atomsFromPDB(pdb2)
contacts2 = getResidueContacts(atoms2, distcut, chainsep)

n=getLength(atoms2)

plotCompContactMap(n, contacts1, contacts2, title = "%s vs %s; cutoff=%s; csep=%i; nc=%i,%i"%(pdb1,pdb2,str(round(distcut,2)),chainsep,len(contacts1),len(contacts2)), fname="%s_VS_%s_contacts.png"%(pdb1[:-4],pdb2[:-4]))

