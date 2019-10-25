import sys
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


def atomsFromGRO(atoms_gro_fname):
	atoms_gro = {}
	gro = open(atoms_gro_fname)
	title = gro.readline()
	natoms = int(gro.readline())
	for line in gro.readlines():
		
		if len(line) >= 44:

			resnr = line[:5].strip()
			residue = line[5:10].strip()
			atype = line[10:15].strip()
			anr = line[15:20].strip()
			x = line[20:28].strip()
			y = line[28:36].strip()
			z = line[36:44].strip()
			atoms_gro[int(anr)] = {"resnr":int(resnr), "residue":residue.strip(), "atype":atype.strip(), "x":float(x), "y":float(y), "z":float(z)}
		
	assert natoms == len(atoms_gro.keys())
	return atoms_gro

def getRes2AtomsDict(atoms):
	resdict = {}
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in resdict:
			resdict[rnum].append(atoms[a])
		else:
			resdict[rnum] = [atoms[a]]
	return resdict

def readXMLContacts(filename):
	contacts = []
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
			ctype = tmp[2]
			parms = map(float,tmp[3:])
			contacts.append([r1,r2, ctype, parms])
	return contacts

def getAtomIDforResidueContact(resid,atype,atoms):
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype:
			return a
	return None


if __name__=="__main__":
	infilename = sys.argv[1]
	contacts = readXMLContacts(infilename)
	outfilename = sys.argv[2]
	pdbfilename = sys.argv[3]
	outfile = open(outfilename, "w")
	atoms = atomsFromPDB(pdbfilename)
	res = getRes2AtomsDict(atoms)

	for c in contacts:
		
		aid1 = getAtomIDforResidueContact(c[0],"CA",atoms) 
		aid2 = getAtomIDforResidueContact(c[1],"CA",atoms)
		print c,aid1,aid2
		if c[2]=="FMULTIGSMOOTH":
			dmin = c[3][0]
			dmax = c[3][1]
			d = (dmin+dmax)/2.0
			w = c[3][4]
			eps = c[3][5]
			outfile.write("\t".join(map(str, [aid1,aid2,eps,d,w]))+"\n")
		elif c[2]=="FMULTIGAUSS":
			d = c[3][0]
			w = c[3][3]
			eps = c[3][4]
			outfile.write("\t".join(map(str, [aid1,aid2,eps,d,w]))+"\n")
	
	outfile.close()