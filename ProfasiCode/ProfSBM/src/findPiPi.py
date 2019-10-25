import sys,math,os,os.path,shutil
import numpy as np
from Bio.PDB.DSSP import *
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning

ringatoms = {"TYR":["CG","CD1","CD2","CE1","CE2","OH"], 
			"PHE":["CG","CD1","CD2","CE1","CE2","CZ"], 
			"TRP":["CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"],
			"HIS":["CG","ND1","CE1","NE2","CD2"],
			"W6C":["CE2","CZ2","CH2","CZ3","CE3","CD2"]}

cations = {"ARG":["CZ","NE"], "LYS":["NZ","CE"]}



def fetchPDB(pdbid, pdb_location,alwaysOverwrite=False,verbose=False):
	success = False
	
	dest = pdb_location + pdbid.upper() + ".pdb"
	
	if not os.path.exists(dest) or alwaysOverwrite:
		if alwaysOverwrite:
			if os.path.exists(dest):
				os.remove(dest) 
		sub = subprocess.Popen("wget -q -P "+pdb_location+" http://www.pdb.org/pdb/files/"+pdbid.upper()+".pdb",shell=True, stdout=subprocess.PIPE)
		#sub2 = subprocess.call("mv "+PDB1.upper()+".pdb "+pdb_location,shell=True, stdout=subprocess.PIPE)
		out = sub.stdout.read()
		#if verbose: print "wget done *****",sub.stdout.read()
		sub.wait()
		
		
	else:
		if verbose:
			print "already found:", pdbid, pdb_location
	
	models = splitModels(dest)
	
	if len(models) == 0:
		print "Could not extract models! ",dest
		return False
	#print models[0], dest
	
	if models[0] != dest:
		shutil.copy(models[0], dest)
	
	if not os.path.exists(dest):
		success = False
	#test basic file integrity
	else:
		for line in open(dest).readlines():
			if line[:4] == "ATOM":
				success = True
				return True
				break
	
	
	return success

def createChainPDBfile(PDBfile, chain, alwaysOverwrite=False):
	pdb1_chain_loc = PDBfile[:-4]+str(chain)+".pdb"
	
	#print "NEW FILE:"+pdb1_chain_loc
	if not alwaysOverwrite and os.path.exists(pdb1_chain_loc):
		#print "EXISTS"
		chain_file = open(pdb1_chain_loc)
		formatOK = True
		linesProcessed = 0
		for line in chain_file.readlines():
			linesProcessed += 1
			if line[:4]!="ATOM":
				formatOK = False
				
				break
				
		chain_file.close()
		if formatOK and linesProcessed > 0:
			return True
	
	
	chain_file = open(pdb1_chain_loc, "w")
	nLinesWritten=0
	for line in open(PDBfile).readlines():
		if line[:4]=="ATOM":# and len(line)>=22:
			if line[21].strip()==chain.strip():
				chain_file.write(line)
				nLinesWritten += 1
			elif line[21].strip()=="":
				modline = line[:21]+"0"+line[22:]
				chain_file.write(modline)
				nLinesWritten += 1
	chain_file.close()
	#print "LINES WRITTEN:",nLinesWritten
	return nLinesWritten > 0

def splitModels(pdbfilename):
	models = []
	model_nr = 0
	atomlines = []
	
	for line in open(pdbfilename).readlines():
		linetype = line[:6]
		
		if linetype == "MODEL ":
			model_nr = int(line.replace("MODEL",""))
			#print "model nr ",model_nr
			atomlines = []
		elif linetype == "ENDMDL":
			assert model_nr != 0, model_nr
			modelfilename = pdbfilename[:-4] + "_m%i.pdb"%model_nr
			model_exists = False
			models.append(modelfilename)
			#if os.path.exists(modelfilename):
			#	model_exists = True
			
			if not model_exists:
				modelfile = open(modelfilename,"w")
				assert len(atomlines) != 0
				for l in atomlines:
					#print l
					modelfile.write(l)
				modelfile.write("END   ")
				modelfile.close()
			
		elif linetype == "ATOM  ":
			atomlines.append(line)
		else:
			pass
	
	if model_nr==0:
		models.append(pdbfilename)
	
	return models # returns a list of the new filenames, one model per file

def check(res,resname):
	ok = True
	
	if resname in ringatoms:
		l = ringatoms[resname]
	elif resname in cations:
		l = cations[resname]
	else:
		print "ERROR: unexpected residues!"
		return False
	
	for a1 in l:
		found = False
		for r in res:
			if r["aname"]==a1:
				found = True
		if not found:
			ok = False
	return ok

def getAtomPoint(res,aname,atoms):
	for at in res:
		if aname==at["aname"]:
			return np.array([at["x"], at["y"], at["z"]])
	return None
	
def atomDist(a1,a2,atoms):
	d = math.sqrt( math.pow(atoms[a1]["x"] - atoms[a2]["x"], 2)
					+ math.pow(atoms[a1]["y"] - atoms[a2]["y"], 2)
					+ math.pow(atoms[a1]["z"] - atoms[a2]["z"], 2)
					)
	#print d
	return d

def dist(a1,a2):
	return math.sqrt( math.pow(a1[0] - a2[0], 2)
					+ math.pow(a1[1] - a2[1], 2)
					+ math.pow(a1[2] - a2[2], 2)
					)
def getAtomIDforResidueContact(resid,atype,atoms):
	
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype:
			
			return a
	
	
	return None

def getResidueDict(atoms):
	rd = {}
	
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in rd:
			assert not atoms[a]["aname"] in [b["aname"] for b in rd[rnum]], "%s\n%s\n%s"%(inpdb,str(rd[rnum]), str(atoms[a]["aname"]))
				
				
			rd[rnum].append(atoms[a])
		else:
			rd[rnum]=[atoms[a]]
	
	
	return rd

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
	prev_rnum = -999999
	for line in open(filename).readlines():
		if line[:4]=="ATOM":
			a = parseATOMline(line)
			rnum=a["rnum"]
			insert = a["insert"]
			altloc = a["altloc"]
			if insert not in [" ", ""]:
				continue
			if altloc not in ["A"," ",""]:
				continue
			if rnum < prev_rnum:
				print inpdb, " - contains multiple conformers!!!"
				break
			atoms[a["anum"]] = a
			prev_rnum = rnum
	return atoms






def getClosestAtomDist(res1,res2,atoms,exclude=""):
	mindist = 10e10
	minat1 = None
	minat2 = None
	for ai in res1:
		if ai["aname"] in exclude: continue
		for aj in res2:
			
			if aj["aname"] in exclude: continue
			d = atomDist(ai["anum"], aj["anum"],atoms)
			if d < mindist:
				mindist = d
				minat1 = ai["anum"]
				minat2 = aj["anum"]
				
	return (mindist,minat1,minat2)
			
	
def getRingCentroid(res,resname):
	#assert resname in ringatoms
	x=[]
	y=[]
	z=[]
	
	for at in res:
		
		name = at["aname"]
		
		if name in ringatoms[resname]:
			#print resname,name,at["x"],at["y"],at["z"]
			x.append(at["x"])
			y.append(at["y"])
			z.append(at["z"])
	#print (np.mean(x), np.mean(y), np.mean(z))
	return np.array([np.mean(x), np.mean(y), np.mean(z)])

def getRingCentroid2(res,resname):
	#assert resname in ringatoms
	x=[]
	y=[]
	z=[]
	points = []
	for at in res:
		
		name = at["aname"]
		
		if name in ringatoms[resname]:
			#print resname,name,at["x"],at["y"],at["z"]
			x = float(at["x"])
			y = float(at["y"])
			z = float(at["z"])
			points.append(np.array([x,y,z]))
	#print (np.mean(x), np.mean(y), np.mean(z))
	points = np.array(points)
	#print "***",points,"***"
	return np.mean(points,axis=0)

def getRingNormal(res):
	
	CG = None
	CD1 = None
	CD2 = None
	for at in res:
		name = at["aname"]
		if name == "CG": CG = np.array([at["x"],at["y"],at["z"]])
		elif name == "CD1": CD1 = np.array([at["x"],at["y"],at["z"]])
		elif name == "CD2": CD2 = np.array([at["x"],at["y"],at["z"]])
	
	if not (CG!=None and CD1!=None and CD2!=None):
		return None
	
	v1 = CD1-CG
	v2 = CD2-CG
	
	#u1 =CG-CD1
	#u2 = CG-CD2
	
	return np.cross(v1,v2)
		

def unit_vector(vector):
	""" Returns the unit vector of the vector.  """
	return vector / np.linalg.norm(vector)

def angle_between_unit(v1, v2):
	""" Returns the angle in radians between vectors 'v1' and 'v2'::
            
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
	"""
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	
	angle = np.arccos(np.dot(v1_u, v2_u))
	if np.isnan(angle):
		if (v1_u == v2_u).all():
			return 0.0
		else:
			return np.pi
	return angle

def angle_between(v1, v2):
	angle = np.arccos(np.dot(v1, v2))
	if np.isnan(angle):
		if (v1 == v2).all():
			return 0.0
		else:
			return np.pi
	return angle 

def r2d(radians):
	if radians != None:
		return radians * (180/np.pi)
	else:
		 return None


def dot_product(x, y):
    return sum([x[i] * y[i] for i in range(len(x))])

def norm(x):
    return math.sqrt(dot_product(x, x))

def normalize(x):
    return [x[i] / norm(x) for i in range(len(x))]

def project_onto_plane(x, n):
    d = dot_product(x, n) / norm(n)
    p = [d * normalize(n)[i] for i in range(len(n))]
    return [x[i] - p[i] for i in range(len(x))]

if __name__=="__main__":
	inpdb = sys.argv[1]
	
	if len(sys.argv)==3:
		givenpath = sys.argv[2]
	else:
		givenpath = None
	
	alwaysOverwrite=False
	
	verbose = False
	
	useW6C = True
	
	pires = ["TRP","TYR","PHE","HIS","ARG","LYS","W6C"]
	
	minClosestDist = 6.0
	
	pdblist = []
	
	if inpdb[-4:] in [".PDB",".pdb"]:
		pdblist.append(inpdb)
	elif givenpath != None:
		assert os.path.isdir(givenpath)
		if os.path.exists(inpdb):
			for line in open(inpdb).readlines():
				tmp = line.strip().split()
				if len(tmp)==1:
					pdblist.append(givenpath+"/"+tmp[0])
	else:
		print "downloading pdb files from list ..."
		c=0
		lines = open(inpdb).readlines()
		for line in lines:
			sys.stdout.write("\r%s%%" % str(round(float(100.0*c)/len(lines),2)))
			sys.stdout.flush()
			c+=1
			
			if len(line.strip()) >= 4:
				code = line.strip()[:4].upper()
				if len(line.strip()) == 4:
					chain = ""
				else:
					chain = line.strip()[-1]
				fname = "./pdb/%s.pdb"%(line.strip())
				if os.path.exists(fname) and not alwaysOverwrite:
					pdblist.append(fname)
				else:
					assert fetchPDB(code,"pdb/",alwaysOverwrite=alwaysOverwrite), ",".join(code,"pdb/",str(alwaysOverwrite))
					assert createChainPDBfile("pdb/"+code+".pdb",chain,alwaysOverwrite=alwaysOverwrite), ",".join("pdb/"+code+".pdb",chain,str(alwaysOverwrite))
					pdblist.append(fname)
			
		
	
	warnings.simplefilter('ignore', BiopythonWarning)
	
	
	print "searching for cation/pi-pi interactions..."
	c = 0
	for inpdb in pdblist:
		
		if not verbose:
			sys.stdout.write("\r%s%%" % str(round(float(100.0*c)/len(pdblist),2)))
			sys.stdout.flush()
			c+=1
		atoms = atomsFromPDB(inpdb)
		residues = getResidueDict(atoms)
		keys = residues.keys()
		if useW6C:
			w6clabel = "_W6C"
		else:
			w6clabel = ""
		outname = inpdb[:-4]+"%s.pidat"%(w6clabel)
		outfile = open(outname,"w")
		length = len(keys)
		
		outfile.write(str(length)+"\n")
		if verbose: print inpdb, length
		#print residues
		try:
			
			structure = PDBParser.get_structure(PDBParser(),inpdb,inpdb)
			model = structure
			dssp_dict, k = dssp_dict_from_pdb_file(inpdb)
		except Warning:
			print "BioPython warned about",inpdb
			#dssp_dict = None
		except Exception:
			print "BioPython doesn't like",inpdb
			dssp_dict = None
		#surface = get_surface(inpdb, PDB_TO_XYZR='./pdb_to_xyzr', MSMS='./msms')

		#rd = ResidueDepth(model, inpdb)

		
		
		datline = "\t".join(["i","j","r1","r2","cdist","mdist","theta","phi","psi","zeta","facc1","facc2"])
		if verbose: print datline
		outfile.write(datline+"\n")
	
		for i in keys:
			ri = residues[i]
			riname = ri[0]["rname"]
			if useW6C and riname=="TRP":
				riname="W6C"
				
			for j in keys:
				rj = residues[j]
				rjname = rj[0]["rname"]
				if useW6C and rjname=="TRP":
					rjname="W6C"
				if i!=j:
					if riname in pires and rjname in pires:
					
						#if riname in cations: continue # and rjname in cations: continue
						
						if dssp_dict != None:
						
							try:
								aa, ss, acc1, phi, psi =  dssp_dict[(chain, (' ',i,' '))]
								aa, ss, acc2, phi, psi =  dssp_dict[(chain, (' ',j,' '))]
								acc1 = acc1/float(MAX_ACC[riname])
								acc2 = acc2/float(MAX_ACC[rjname])
							except Exception:
								acc1 = None
								acc2 = None
						else:
							acc1 = None
							acc2 = None
						
						min_dist = getClosestAtomDist(ri,rj,atoms,exclude=["CA","C","O","N","CB","H","H1","H2","H3","HA","HB","HB1","HB2","HB3","HD","HD1","HD2","HD3","HE","HE1","HE2","HE3","HG","HG1","HG2","HG3","HH","HZ","HZ1","HZ2","HZ3"])[0]
						
						
						
						if riname in ringatoms:
							centroid_i = getRingCentroid(ri,riname)
							
								
							#print centroid_i, getRingCentroid2(ri,riname)
							#assert centroid_i.all() == getRingCentroid2(ri,riname).all()
						else:
							centroid_i = getAtomPoint(ri,cations[riname][0],atoms)
						
						#if centroid_i == None:
						#	print "Error:",i,riname,cations[riname][0],[x["aname"]for x in ri]
						#	continue
						
						if rjname in ringatoms:
							centroid_j = getRingCentroid(rj,rjname)
						else:
							centroid_j = getAtomPoint(rj,cations[rjname][0],atoms)
						
						if not check(ri,riname): centroid_i = None
						if not check(rj,rjname): centroid_j = None
						
						if centroid_j == None or centroid_i == None:
							#print "Error:",j,rjname,cations[rjname][0],[x["aname"]for x in rj]
							datline = "\t".join([str(round(d,3)) if type(d) in [float,np.float64] else str(d) for d in [i,j,riname,rjname,None,min_dist,None,None,None,None,acc1,acc2]])
							if verbose: print datline
							outfile.write(datline+"\n")
						else:
					
							cent_dist = np.linalg.norm(centroid_i-centroid_j)
							#cent_dist2 = math.sqrt( math.pow(centroid_j[0]-centroid_i[0],2) + math.pow(centroid_j[1]-centroid_i[1],2) + math.pow(centroid_j[2]-centroid_i[2],2) )
							#print cent_dist,cent_dist2
							#print cent_dist, np.linalg.norm(unit_vector(centroid_i) - unit_vector(centroid_j))
							#assert cent_dist == cent_dist2
							#assert np.linalg.norm(centroid_i-centroid_j) == np.linalg.norm(centroid_j-centroid_i) == np.linalg.norm(unit_vector(centroid_i) - unit_vector(centroid_j))
					
							if riname in ringatoms:
								normal_i = getRingNormal(ri)
							else:
								normal_i = np.array(getAtomPoint(ri,cations[riname][1],atoms) - getAtomPoint(ri,cations[riname][0],atoms)) # instead of normal, use two atoms from cation
					
							if rjname in ringatoms:
								normal_j = getRingNormal(rj)
							else:
								normal_j = np.array(getAtomPoint(rj,cations[rjname][1],atoms) - getAtomPoint(rj,cations[rjname][0],atoms)) # instead of normal, use two atoms from cation
					
							if normal_i != None and normal_j != None:
								theta = angle_between_unit(normal_i, normal_j)
								#print r2d(theta), r2d(np.pi-theta)
								theta = min(theta, np.pi-theta)
							else:
								theta = None
					
							if normal_i != None and normal_j != None:
								phi1 = angle_between_unit(normal_i, centroid_j-centroid_i)
								phi1 = min(phi1, np.pi-phi1)
							else:
								phi1 = None
							
							
					
					
					
							# psi angle: clockwise difference between centroid-CG vectors for two aromatics
						
						
							if riname in ringatoms and rjname in ringatoms:
								xi = getAtomPoint(rj,ringatoms[rjname][0],atoms)
								x = np.array(xi - centroid_j) # CG-centroid vector
								n = normal_i
								v = np.array(getAtomPoint(ri,ringatoms[riname][0],atoms) - centroid_i) # CG-centroid vector
							elif riname in ringatoms and rjname in cations:
								x = np.array(getAtomPoint(rj,cations[rjname][1],atoms) - getAtomPoint(rj,cations[rjname][0],atoms))
								n = normal_i
								v = np.array(getAtomPoint(ri,ringatoms[riname][0],atoms) - centroid_i) # CG-centroid vector
							elif riname in cations and rjname in ringatoms:
								x = np.array(getAtomPoint(rj,ringatoms[rjname][0],atoms) - centroid_j) # CG-centroid vector
								n = None
								v = np.array(getAtomPoint(ri,cations[riname][1],atoms) - getAtomPoint(ri,cations[riname][0],atoms))
							
						
							if n==None:
								psi1 = None
							else:
								proj_j2i = project_onto_plane(x,n)
								psi1 = angle_between_unit(v, proj_j2i)
								#psi1 = min(psi1, 2*np.pi -psi1)
						
						
							if riname in ringatoms and rjname in ringatoms:
								x = centroid_j
								n = normal_i
								v = np.array(getAtomPoint(ri,ringatoms[riname][0],atoms) - centroid_i) # CG-centroid vector
							elif riname in ringatoms and rjname in cations:
								x = centroid_j
								n = normal_i
								v = np.array(getAtomPoint(ri,ringatoms[riname][0],atoms) - centroid_i) # CG-centroid vector
							elif riname in cations and rjname in ringatoms:
								x = centroid_j
								n = None
								v = np.array(getAtomPoint(ri,cations[riname][1],atoms) - getAtomPoint(ri,cations[riname][0],atoms))
						
							if n == None:
								zeta = None
							else:
								proj_j2i = np.array(project_onto_plane(x,n))
								zeta = angle_between_unit(v, proj_j2i-centroid_i)
						
						
						
					
						
						#print MAX_ACC[riname], MAX_ACC[rjname]
					
						#rd = residue_depth(i+1, surface)
						#print i+1,aa,ss,acc#,rd
						if min_dist <= minClosestDist:
							datline = "\t".join([str(round(d,3)) if type(d) in [float,np.float64] else str(d) for d in [i,j,riname,rjname,cent_dist,min_dist,r2d(theta),r2d(phi1),r2d(psi1),r2d(zeta),acc1,acc2]])
							if verbose: print datline
							outfile.write(datline+"\n")
		outfile.close()