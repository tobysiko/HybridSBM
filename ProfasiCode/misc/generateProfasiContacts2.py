import sys
import numpy
import optparse
import math
import subprocess
import os, os.path
import copy

BB_atoms = ["C","CA","O","N","OXT","CB"]

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
					'VAL':'V'
}
AA_1_to_3 = {
					'G':'GLY',
					'A':'ALA',
					'V':'VAL',
					'L':'LEU',
					'I':'ILE',
					'C':'CYS',
					'M':'MET',
					'F':'PHE',
					'Y':'TYR',
					'W':'TRP',
					'P':'PRO',
					'S':'SER',
					'T':'THR',
					'N':'ASN',
					'Q':'GLN',
					'D':'ASP',
					'E':'GLU',
					'H':'HIS',
					'K':'LYS',
					'R':'ARG'
}

aromatic = ['PHE','TYR','TRP']


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

def translateAtoms(a1,a2,trans_atoms,ref_atoms):
	#assert a1 in trans_atoms and a2 in trans_atoms, (a1,a2)
	#assert a1 in ref_atoms and a2 in ref_atoms, (a1,a2)
	trans = []
	# print a1,a2
	for a in [a1,a2]:
		trans_aname = trans_atoms[a]["aname"]
		trans_rname = trans_atoms[a]["rname"]
		trans_rnum = trans_atoms[a]["rnum"]
		trans_chain = trans_atoms[a]["chain"]
		found = False
		for ra in ref_atoms.keys():
			
			ref_aname = ref_atoms[ra]["aname"]
			ref_rname = ref_atoms[ra]["rname"]
			ref_rnum = ref_atoms[ra]["rnum"]
			ref_chain = ref_atoms[ra]["chain"]
			
			if ref_aname==trans_aname and ref_rnum==trans_rnum and ref_chain==trans_chain:
				
				found=True
				#print "found atom",a
				break
		if found:
			trans.append(ra)
		else:
			trans.append(None)
	
	return trans[0], trans[1]
	
	

def filterSuperEnsembleContacts(models_contact_lists, models_atoms, atypes, ensembleCutoff, contactType, verbose=False):
	print "super ensemble contacts"
	contactsDict = {}
	
	atomtypes = copy.copy(atypes)
	
	if contactType == "HYBRID" and not "CA" in atomtypes:
		atomtypes  += ['CA']
	
	nModels = len(models_contact_lists)
	contactsPerModel = [[] for i in xrange(nModels)]
	
	
	
	for m in xrange(nModels):
		model_contacts = models_contact_lists[m]
		print "%i contacts in model %s"%(len(model_contacts),m)
		
		modelatoms = models_atoms[m]
		print len(modelatoms)," atoms in modelatoms"
		
		for c in model_contacts:
			a1 = c[0]
			a2 = c[1]
			d  = c[2]
			
			
			
			
			assert modelatoms[a1]["aname"] in atypes, modelatoms[a1]
			assert modelatoms[a2]["aname"] in atypes, modelatoms[a2]
			#contactsPerModel[m].append( [a1,a2,d,modelatoms] )
			
			
			#ref_a1,ref_a2 = translateAtoms(a1,a2,modelatoms,ref_atoms)
			
			atype1 = modelatoms[a1]["aname"]
			atype2 = modelatoms[a2]["aname"]
			rnum1 = modelatoms[a1]["rnum"]
			rnum2 = modelatoms[a2]["rnum"]
			
			#print ref_atoms[ref_a1], ref_atoms[ref_a2]
			
			
			
			if (rnum1,rnum2,atype1,atype2) in contactsDict:
				contactsDict[(rnum1,rnum2,atype1,atype2)].append(d)
			else:
				contactsDict[(rnum1,rnum2,atype1,atype2)] = [d]
	
	# which contacts are found in all models? what are the means and standard deviations of distances?

	counter = 0
	#avgStdDistsShared = []
	
	finalContacts = []
	dists = []
	#means = []
	#stddevs = []
	
	global_min_dist = 10000
	global_max_dist = -1
	min_range = 100000
	max_range = -1
	
	
	
	print "Collected %i contacts from all models."%len(contactsDict.keys())
	
	for i in sorted(contactsDict.keys()):
		
		contactCount = len(contactsDict[i]) 
		assert contactCount <= nModels, str([contactCount, nModels])
		
		#print contactsDict[i]
		
		frac = contactCount / float(nModels)
		
		if frac >= ensembleCutoff: # in how many models was this contact found ? Is the fraction above threshold?
			counter += 1
			
			#avgDist = numpy.mean(contactsDict[i])
			#stdDist = numpy.std(contactsDict[i])
			
			
			min_d = min(contactsDict[i])
			max_d = max(contactsDict[i])
			range_d = max_d - min_d
			
			if min_d < global_min_dist: global_min_dist = min_d
			if max_d > global_max_dist: global_max_dist = max_d
			if range_d < min_range: min_range = range_d
			if range_d > max_range: max_range = range_d
			
			
			print "(%s) accepted atom pair:"%str(round(frac,1)), i, " - ", contactCount, "models with this contact.", "Min,Max:",min_d,max_d," Range:",range_d
			
			
			
			if i[2] in atomtypes and i[3] in atomtypes:
				#assert atoms[i[0]]["aname"] in atomtypes
				#assert atoms[i[1]]["aname"] in atomtypes
				#print contactsDict[i]
				fc = list(i) + [contactsDict[i]]
				#print fc
				finalContacts.append(fc)
				
				
				#means.append(avgDist)
				#stddevs.append(stdDist)
			
				#avgStdDistsShared.append([i[0], i[1], avgDist, stdDist])
			else:
				print "atom is not allowed for contact!",ref_atoms[i[0]]["aname"],ref_atoms[i[1]]["aname"]
		else:
			print "(%s) REJECTED atom pair:"%str(round(frac,1)), i, " - ", contactCount, "models with this contact."
	
	print "Min and Max distances within all contacts:", global_min_dist, global_max_dist
	print "Min and Max distance ranges:", min_range,max_range
	
	print "%i contacts above threshold"%(counter)
	
	
	
	return finalContacts#, means, stddevs




def getCutoffContacts(atoms, atomtypes, cutoff, chainDist, verbose=False, res2SS = None, selection=[], selectExclusive=False):
	#atoms = atomsFromPDB(pdbfilename)
	#print atomtypes
	
	#print len(atoms), "atoms"
	
	contactAtoms = []
	
	for a in sorted(atoms.keys()):
		#print a, atoms[a]
		if atoms[ a ][ "aname" ] in atomtypes:
			contactAtoms.append( a )
	
	contacts = []

	for i in xrange( len( contactAtoms ) ):
		ai = contactAtoms[i]
		ri = atoms[ ai ][ "rnum" ]
		if selectExclusive and not ri in selection: continue
		if res2SS != None and res2SS[ri][0] == 'X':
			#print "skip ri",ri,res2SS[ri]
			continue

		for j in xrange( len( contactAtoms ) ):
			aj = contactAtoms[j]
			rj = atoms[ aj ][ "rnum" ]
			if selection!=[] and selectExclusive and not rj in selection: continue
			
			if selection!=[] and not selectExclusive and not ri in selection and not rj in selection: continue
			
			if res2SS != None and res2SS[rj][0] == 'X':
				#print "skip rj",rj,res2SS[rj]
				continue
			
			if ai < aj:
				if ri < rj - chainDist:
					d = atomDist( ai, aj, atoms )
					#print d
					if d <= cutoff:
						contacts.append( [ai, aj, d] )
						if verbose: print [ai, aj, d], ri,rj
	#print contacts
	
	print "Found %i CUTOFF contacts between %s atoms at chain separation=%i and cutoff distance = %f"%(len(contacts), str(atomtypes),chainDist, cutoff)
	
	return contacts

def getResidueDict(atoms):
	rd = {}
	
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in rd:
			rd[rnum].append(atoms[a])
		else:
			rd[rnum]=[atoms[a]]
	
	
	return rd

def getResidueContacts(atoms, cutoff, chainDist, verbose=False, res2SS = None, selection=[], selectExclusive=False):
	#from Bio.PDB import *
	
	#structure = PDBParser.get_structure(PDBParser(),pdbfilename,pdbfilename)

	cutoff = cutoff

	chainsep = chainDist

	#model = structure[0]

	#reslist = Selection.unfold_entities(structure, 'R')

	contacts = []
	residues = getResidueDict(atoms)
	n=len(residues)
	print n,"residues"
	
	for i in xrange(n):
		if selectExclusive and not i+1 in selection: continue
		ri = residues[i+1]
		for j in xrange(n):
			if selection!=[] and selectExclusive and not j+1 in selection: continue
			
			if selection!=[] and not selectExclusive and not i+1 in selection and not j+1 in selection: continue
			
			rj = residues[j+1]
		
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
					
					contacts.append( [a1,a2,cadist] )
					#print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	
	
	return contacts
def getAtomIDforResidueContact(resid,atype,atoms):
	
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype:
			
			return a
	
	
	return None

def gro2pdbAtomMap(grofilename, pdbfilename, atomtypes):
	gro2pdb_heavy = {}
	groatoms = atomsFromGRO(grofilename)
	pdbatoms = atomsFromPDB(pdbfilename)
	
	#print len(groatoms), len(pdbatoms)
	
	assert (len(groatoms) - len(pdbatoms)) <= 1
	
	for i in groatoms.keys():
		resnr = groatoms[i]["resnr"]
		atype = groatoms[i]["atype"]
		if not atype in atomtypes: continue
		found = False
		for j in pdbatoms.keys():
			rnum = pdbatoms[j]["rnum"]
			aname = pdbatoms[j]["aname"]
			if atype == "OC1" and aname == "OXT":
				atype = "OXT"
			
			
			if resnr == rnum and atype == aname:
				gro2pdb_heavy[i] = j
				found = True
		if not found:
			print "could not find",i,resnr,atype
	
	return gro2pdb_heavy
	


def getShadowContacts(pdbfilename, outfilename, cutoff=6.0, chainDist=3, gmxpath="", scmpath="", res2SS = None, verbose = False, selection=[], selectExclusive=False):
	
	cleanupList = []
	
	grofilename = pdbfilename[:-4]+".gro"
	topfilename = pdbfilename[:-4]+".top"
	
	pdb2gmx = "%spdb2gmx -f %s -o %s -p %s -ff amber99sb-ildn -water tip3p -ignh -i tmpposre.itp"%(gmxpath, pdbfilename, grofilename, topfilename)
	if verbose: print pdb2gmx
	sub = subprocess.Popen(pdb2gmx,stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	
	#print err
	
	atoms = atomsFromGRO( grofilename )
	
	cleanupList.append(grofilename)
	cleanupList.append(topfilename)
	cleanupList.append("tmpposre.itp")
	
	
	scm = "java -jar %sSCM.jar -t %s -g %s -o %s --distance --coarse AACA -m shadow -s 1.0 -c %s --proteinDelta %s"%(scmpath, topfilename, grofilename, outfilename, str(cutoff), str(chainDist))
	if verbose: print scm
	sub = subprocess.Popen(scm,stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	contacts = []
	
	#print err
	
	gro2pdb = gro2pdbAtomMap(grofilename, pdbfilename, atomtypes=['CA'])
	
	for line in open(outfilename).readlines():
		tmp = line.strip().split()
		if len(tmp)==5:
			a1 = int(tmp[1])
			a2 = int(tmp[3])
			d = float(tmp[4])*10.0
			r1 = atoms[ a1 ][ "resnr" ]
			r2 = atoms[ a2 ][ "resnr" ]

			if selection!=[] and not selectExclusive and (not r1 in selection or not r2 in selection): continue
			
			if selection!=[] and selectExclusive and not r1 in selection and not r2 in selection: continue
			#print atoms[ a1 ]
			#print atoms[ a2 ]
			
			a1 = gro2pdb[a1]
			a2 = gro2pdb[a2]
			if res2SS != None:
				if res2SS[r1][0] != 'X' and res2SS[r2][0] != 'X':
					contacts.append( [ a1, a2, d ] )
			else:
				contacts.append( [ a1, a2, d ] )
	
	cleanupList.append(outfilename)
	
	for i in cleanupList:
		if os.path.exists(i):
			os.remove(i)
	
	print "Found %i SHADOW contacts in %s between CA atoms at chain separation=%i and cutoff distance = %f"%(len(contacts), pdbfilename, chainDist, cutoff)
	return contacts


def atomDist(a1,a2,atoms):
	d = math.sqrt( math.pow(atoms[a1]["x"] - atoms[a2]["x"], 2)
					+ math.pow(atoms[a1]["y"] - atoms[a2]["y"], 2)
					+ math.pow(atoms[a1]["z"] - atoms[a2]["z"], 2)
					)
	#print d
	return d



def mixCutoffShadow(contacts_co, contacts_sh):
	contacts = []
	
	for c in contacts_co: contacts.append(c)

	for sh in xrange(len(contacts_sh)):
		a1_sh, a2_sh, d_sh = contacts_sh[sh]
		found = False
		for co in xrange(len(contacts_co)):
			a1_co, a2_co, d_co = contacts_co[co]
			
			if a1_sh == a1_co and a2_sh == a2_co: # shadow contact is not in cutoff contacts
				found = True
				break
		if not found: 
			contacts.append(contacts_sh[sh])
	return contacts


def profasiAtomString(anum, atoms):
	rnum = atoms[anum]["rnum"]
	aname = atoms[anum]["aname"]	
	rname = atoms[anum]["rname"]
	return "0/%s/%s/_%s_"%(str(rnum-1), rname, aname)

def profasiAtomStringNoResType(anum, atoms):
	rnum = atoms[anum]["rnum"]
	aname = atoms[anum]["aname"]	
	rname = atoms[anum]["rname"]
	return "0/%s//_%s_"%(str(rnum-1),  aname)


def atoms2atoms(atoms1, atoms2, atomtypes):
	a2a = {}
	for i in atoms1.keys():
		resnr = atoms1[i]["rnum"]
		atype = atoms1[i]["aname"]
		if not atype in atomtypes: continue
		found = False
		for j in atoms2.keys():
			rnum = atoms2[j]["rnum"]
			aname = atoms2[j]["aname"]
			
			if resnr == rnum and atype == aname:
				a2a[i] = j
				found = True
		if not found:
			print "could not find",i,resnr,atype
	return a2a

def getResNameFromNumber(atoms, rnum, aname):
	for a in atoms:
		if atoms[a]["rnum"]==rnum and atoms[a]["aname"]:
			return atoms[a]["rname"]
	
	return None

def getRes2AtomsDict(atoms):
	resdict = {}
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in resdict:
			resdict[rnum].append(atoms[a])
		else:
			resdict[rnum] = [atoms[a]]
	return resdict

def writeProfasiContactsFMULTIGAUSS(contacts, atoms, atomtypes, outfilename, radius, steepness, width, depth, norm=0,label="",includeNonNative=False, chainSep=3, smoothGaussian=False, restraintType="distance_restraints"):
	
	
	
	n = len(contacts)
	if norm > 0:
		depth = norm/float(n)
		print "New well depth:",depth
	
	outfilename = outfilename + "_%icontacts_E%s"%(n, str(round(n*depth,1))) + ".xml"
	
	out = open(outfilename,"w")
	
	if smoothGaussian:
		out.write("""<%s rid="%s">  
<formatted_data>"
    <format name=\"restraint\" type=\"$3\">
      <atom1>$1</atom1>"
      <atom2>$2</atom2>
      <parameters>
        <low>$4</low>
        <high>$5</high>
        <radius>$6</radius>
        <steepness>$7</steepness>
        <width>$8</width>
        <depth>$9</depth>
      </parameters>
    </format>
<data>"""%(restraintType,label))
	else:
		out.write("""<%s rid="%s">  
<formatted_data>"
    <format name=\"restraint\" type=\"$3\">
      <atom1>$1</atom1>"
      <atom2>$2</atom2>
      <parameters>
        <minima>$4</minima>
        <radius>$5</radius>
        <steepness>$6</steepness>
        <width>$7</width>
        <depth>$8</depth>
      </parameters>
    </format>
<data>"""%(restraintType,label))
	
	seqdists = 0
	#print radius, steepness, width, depth

	for ci in xrange(len(contacts)):
		c = contacts[ci]
		rnum1 = c[0]
		rnum2 = c[1]
		seqdists += abs(rnum1-rnum2)
		#print rnum1,rnum2
		atype1 = c[2]
		atype2 = c[3]
		dists = c[4]
		energies = [getEnergy(i,dists,radius, steepness, width, depth) for i in dists]
		assert all([i==-depth for i in energies])
		
		mindist = min(dists)
		if radius < 0:
			final_radius = mindist + radius
		else:
			final_radius = radius
		#print rnum1,rnum2,atype1,atype2
		
		
		
		atomstring1 = "0/%s/%s/_%s_"%(str(rnum1-1), getResNameFromNumber(atoms,rnum1,atype1), atype1)#profasiAtomString(a1, atoms)
		atomstring2 = "0/%s/%s/_%s_"%(str(rnum2-1), getResNameFromNumber(atoms,rnum2,atype2), atype2)#profasiAtomString(a2, atoms)
		
		if smoothGaussian:
			out.write( "  %s  %s  FMULTIGSMOOTH  %s  %s  %s %s %s %s\n"%(atomstring1, atomstring2, str(min(dists)), str(max(dists)), str(final_radius), str(steepness), str(width), str(depth) ) )
		else:	
			out.write( "  %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n"%(atomstring1, atomstring2, ",".join([str(i) for i in dists]), str(final_radius), str(steepness), str(width), str(depth) ) )
		
	
	N=getLength(atoms)
	print "CO = ", seqdists/float(n*N), seqdists, n, N
	
	
	if includeNonNative:
		nnCounter = 0
		nCounter = 0
		resdict = getRes2AtomsDict(atoms)
		#N = len(resdict)
		
		#print "CO = ", seqdists/float(n*N)
		
		for i in xrange(N):
			ai = None
			for a in resdict[i+1]: 
				if a["aname"] in atomtypes: ai = a
			assert ai != None
			
			for j in xrange(N):
				assert isinstance(i ,(int))
				if i < j - chainSep:
					aj = None
					for a in resdict[j+1]: 
						if a["aname"] in atomtypes: aj = a
					assert aj != None
					
					isNative = False
					
					
					#print (i,j)
					for ci in xrange(len(contacts)):
						c = contacts[ci]
						rnum1 = c[0]
						rnum2 = c[1]
						atype1 = c[2]
						atype2 = c[3]
						if rnum1==i+1 and rnum2==j+1 and ai["aname"]==atype1 and aj["aname"]==atype2:
							#print i+1,j+1,ai["aname"],aj["aname"],"yes"
							isNative = True
						#else:
						#	#print i+1,j+1,ai["aname"],aj["aname"],"no"
					
					if not isNative:
						nnCounter += 1
						atomstring1 = "0/%s/%s/_%s_"%(str(i), ai["rname"], ai["aname"])
						atomstring2 = "0/%s/%s/_%s_"%(str(j), aj["rname"], aj["aname"])
		
						out.write( "  %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n"%(atomstring1, atomstring2, ",".join(["0" for x in dists]), "4", "1", "0", "0" ) )
					else:
						nCounter += 1
		print nnCounter, "non-native contacts written! (chain sep = %i)"%chainSep
		print nCounter, "native contacts"
		assert nCounter == n
		assert nnCounter+nCounter == (N*N)/2.0 - N*0.5 - sum([N-xx for xx in range(1,chainSep+1)]), str(nnCounter+nCounter)+ " "+str( (N*N)/2.0 - N*0.5 - sum([N-xx for xx in range(1,chainSep+1)])  )
	out.write("""    </data>
  </formatted_data>
</%s>
"""%(restraintType))
	print "Expected energy:",n*depth
	
	return outfilename


def writeChimeraPseudobonds(contacts, outfilename, legend=None):

	outfile = open(outfilename, 'w')
	
	if legend != None: assert len(contacts) == len(legend),str(len(contacts))+" "+str(len(legend))
	
	maxrange = 0
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		dists = c[4]
		rn = max(dists) - min(dists)
		if rn > maxrange:
			maxrange = rn
	print "maxrange",maxrange
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		rnum1 = c[0]
		rnum2 = c[1]
		atype1 = c[2]
		atype2 = c[3]
		dists = c[4]
		rn = float(max(dists) - min(dists))
		
		label = str(round(min(dists),2)) + "-" + str(round(max(dists),2))
		
		if legend==None:
			#color = "blue"
			if maxrange > 0:
				gray = int(math.ceil((1.0-(rn/maxrange)) * 65535))
			else:
				gray = 65535
			#print gray
			hx = hex(gray)
			color = "#"+str(hx).replace("0x","")*3
			#color = color.upper()
			#print gray,color,hx
		else:
			lstring = legend[ci]
			if lstring == "s1":
				color = "blue"
			elif lstring == "s2":
				color = "red"
			else:
				color = "magenta"
		pseudobond = "#0:%s@%s #0:%s@%s %s %s"%(str(rnum1), atype1.lower(), str(rnum2), atype2.lower(), color, label)
		outfile.write(pseudobond + "\n")

	outfile.close()




def getEnergy(x, minima, radius, steepness, width, depth):
	
	wells = 1.0
	for m in minima:
		if width != 0.0:
			Gij = - math.exp( - math.pow(x-m, 2.0) / (2 * width * width) )
	
		else:
			Gij = 0.0
	
		wells *=  (1 + Gij)

	Rij = steepness * math.pow(radius/x, 12.0)
	repulsion = (1 + (Rij/depth) )
		
	return (depth * repulsion * wells ) - depth


def plotContactMap(n, contacts, inPDBFilename, cutoff, chainSep, atoms, contactType, fname="contactmap.png"):
	
	import matplotlib.pyplot as p
	import matplotlib.cm as cm
	
	cdict = {}
	for c in contacts:
		r1 = c[0]
		r2 = c[1]
		dists = c[4]
		
		#cdict[(r1,r2)] = numpy.mean(dists) # plot MEAN CA distances
		if len(dists)==1:
			cdict[(r1,r2)] = 1.0
		else:
			cdict[(r1,r2)] = max(dists)-min(dists)  # plot RANGE of CA distances
	
	matrix = [["NaN" for j in xrange(n+1)] for i in xrange(n+1) ]
	
	
	
	#cmap = cm.binary
	for i in xrange(n+1):
		
		for j in xrange(n+1):
			
			
			
			if (i,j) in cdict:
				matrix[i][j] = cdict[(i,j)]
				matrix[j][i] = cdict[(i,j)]
				#p.scatter([i+1],[j+1],marker="s",s=5,cmap=cmap)
			#else:
			#	matrix[i][j] = None
			#	matrix[j][i] = None
		
	
	open(fname+".data","w").write("\n".join( [",".join([str(j) for j in i[1:]]) for i in matrix[1:]]))
	
	matrix = [[matrix[i][j] if matrix[i][j] != "NaN" else 0.0 for j in xrange(n+1)] for i in xrange(n+1) ]
	
	title = "%s; cutoff=%s; csep=%i; n=%i\n%s"%(inPDBFilename,str(round(cutoff,2)),chainSep,len(contacts),contactType)

	imax = p.matshow(matrix,cmap=cm.afmhot,aspect='equal', origin='lower')
	
	p.plot([1,n],[1,n],color='k')
	p.title(title)
	cbar = p.colorbar()
	p.xlabel("residue position")
	p.ylabel("residue position")
	p.xlim(0.5,n+0.5)
	p.ylim(0.5,n+0.5)
	p.grid(color="0.5", linestyle=':', linewidth=1.0)
	#cbar.set_label("CA dist")
	cbar.set_label("CA dist. range")
	#imax.get_axes().set_axisbelow(True)
	p.savefig(fname,dpi=600)
	
	#p.show()
	print "wrote contact map PNG:",fname
	






def getLength(atoms):
	maxrnum = 0
	for a in atoms:
		if atoms[a]["rnum"]>maxrnum:
			maxrnum = atoms[a]["rnum"]
	return maxrnum



	

def getPrfSCcontactEnergy(tconf, pdb, selection, parameters, prf_loc="~/workspace/PROFASI/app/bin/"):
	#~/workspace/PROFASI/app/bin/prf_energies 1PGA_min_rmsd.tconf --box_length 1000 --skip_energy FF08:ExVol,FF08:Bias,FF08:LocExVol,FF08:HBMM,FF08:TorsionTerm --add_chain_pdb 1 1PGA_min_rmsd.pdb::A,45,45:A,49,49 -force_field FF14 --cationpi 10
	cmd = "%sprf_energies %s --box_length 1000 --skip_energy FF08:ExVol,FF08:Bias,FF08:LocExVol,FF08:HBMM,FF08:TorsionTerm --add_chain_pdb 1 %s%s -force_field FF14 %s" % (prf_loc, tconf, pdb, selection, parameters)
	#print cmd
	sub = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	print err
	print out
	HBMS = None
	Hydrophobicity = None
	ChargedSCInteraction = None
	CationPi = None
	Total = None
	start = False
	for line in out.split("\n"):
		if ">>>>>>>>>>" in line:
			#print "start"
			start = True
			continue
		
		if start:
			tmp = line.strip().split()
			if len(tmp)==3:
				#print tmp
				if tmp[0]=="HBMS": HBMS = float(tmp[2])
				elif tmp[0]=="Hydrophobicity": Hydrophobicity = float(tmp[2])
				elif tmp[0]=="ChargedSCInteraction": ChargedSCInteraction = float(tmp[2])
				elif tmp[0]=="CationPi": CationPi = float(tmp[2])
				elif tmp[0]=="Total": Total = float(tmp[2])
	return (HBMS,Hydrophobicity,ChargedSCInteraction,CationPi,Total)

def getSCenergies4Contacts(contacts,pdb,atoms):
	energies=[]
	print "Profasi SC energies per contact:"
	for c in contacts:
		
		r1 = c[0]
		r2 = c[1]
		a1 = c[2]
		a2 = c[3]
		selection = "::A,%i,%i:A,%i,%i"%(r1,r1,r2,r2)
		parameters = "--cationpi 10"
		tconf = pdb[:-4]+".tconf"
		sc_energy = getPrfSCcontactEnergy(tconf, pdb, selection, parameters, prf_loc="~/workspace/PROFASI/app/bin/")
		print c, sc_energy,getResNameFromNumber(atoms,r1,a1),getResNameFromNumber(atoms,r2,a2)
		energies.append(sc_energy)
	return energies


def multiGaussSmoothFunc(x, dists, eps, w):
	vals = []
	low = min(dists)
	high = max(dists)
	for i in x:
		if i <= low:
			vals.append(eps*( 1.0 - numpy.exp( - (numpy.power(i-low,2.0)/(2.0*w*w)) )) - eps )
		elif i >= high:
			vals.append(eps*( 1.0 - numpy.exp( - (numpy.power(i-high,2.0)/(2.0*w*w)) )) - eps )
		else:
			vals.append( -eps )
	return numpy.array(vals)
def multiGaussFunc(x, dists, eps, w):
	prod = eps
	for i in dists:
		prod *= 1.0 - numpy.exp( - (numpy.power(x-i,2.0)/(2.0*w*w)) ) 
	return eps*prod - eps
def LennardJones(x, dist, eps):
	return eps * (numpy.power(dist/x,12) - 2*numpy.power(dist/x,6) )
def multiLJ(x, dists, eps):
	return sum([LennardJones(x,d,eps) for d in dists])
def plotMultiGaussianContacts(contacts, eps, w, fname, smoothGaussian=False):
	import matplotlib.pyplot as p
	l = len(contacts)
	fig = p.Figure(figsize=(25, 25),dpi=800)
	p.suptitle("")
	sq = int(math.ceil(math.sqrt(len(contacts))))
	
	globmin = 10e10
	globmax = -10e10
	for ci in xrange(len(contacts)):
		c=contacts[ci]
		r1 = c[0]
		r2 = c[1]
		dists = c[4]
		if min(dists) < globmin: globmin = min(dists)
		if max(dists) > globmax: globmax = max(dists)
	
	for ci in xrange(len(contacts)):
		c=contacts[ci]
		r1 = c[0]
		r2 = c[1]
		dists = c[4]
		
		ax = p.subplot(sq, sq, ci+1)
		#ax.set_title("%i: (%i,%i)"%(ci+1,r1,r2),size=4)
		
		x = numpy.arange(math.floor(globmin)-1.0,math.ceil(globmax)+1.0,0.01)
		ax.set_ylim(-eps-0.1,0.2)
		ax.set_xlim(math.floor(globmin)-1.0,math.ceil(globmax)+1.0)
		
		textsize = 5
		
		
		if ci in [i for i in xrange(len(contacts)) if i==0 or i%sq==0 ]: 
			#ax.set_yticklabels([0.2,0.0,-0.2,-0.4,-0.6,-0.8,-1.0],size=textsize)
			#ax.set_yticklabels([i.get_text() for i in ax.get_yticklabels()],size=textsize)
			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(textsize-1)
			ax.set_ylabel("epsilon",size=textsize)
		else: 
			#ax.set_yticks([])
			ax.set_yticklabels([])
		
		if ci in [i for i in range(len(contacts)-sq, len(contacts) ) ]:	
			#ax.set_xticklabels(range(int(math.floor(globmin)-1.0),int(math.ceil(globmax)+1.0),1), size=textsize)
			#ax.set_xticklabels([i.get_text() for i in ax.get_xticklabels()],size=textsize)
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(textsize-1)
			ax.set_xlabel("dij",size=textsize)
		else: 
			#ax.set_xticks([])
			ax.set_xticklabels([])
		
		for d in dists:
			p.plot(x, LennardJones(x, d, eps),linestyle='-',linewidth=0.1,c='red')
		p.plot(x, multiLJ(x, dists, eps)/len(dists),linestyle='-',linewidth=0.4,c='green')
		
		if smoothGaussian:
			p.plot(x, multiGaussSmoothFunc(x, dists, eps, w),linewidth=0.8,c='blue')
		else:
			p.plot(x, multiGaussFunc(x, dists, eps, w),linewidth=0.8,c='blue')
		
		
		p.text(globmin, 0.3, "%i: (%i,%i)"%(ci+1,r1,r2),size=textsize)
	
	p.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.4)
	p.savefig(fname,dpi=600)
	#p.show()

if __name__=="__main__":
	##################
	# HANDLE OPTIONS (-h for help)
	description="%prog prepares native contact files in XML format for PROFASI as modified by TS.\n Default values for options are given in parentheses."
	version="2014-09-17; written by Tobias Sikosek"
	usage="\n\tpython %prog [OPTIONS] -s filename.pdb"
	
	op=optparse.OptionParser(description=description, usage=usage, version=version)
	
	op.add_option('-s','--inPDBFilename', action="store", type="string", dest="inPDBFilename",help="structure file in PDB format. Comma-separate multiple files in order to use potential type FMULTIGAUSS and omit -S option.")
	op.add_option('-c','--contactType', action="store", type="string", dest="contactType",help="type of contact definition: CUTOFF, SHADOW, HYBRID, RESIDUE  (%default)")
	op.add_option('-a','--contactAtoms', action="store", type="string", dest="contactAtoms",help="list of atom types for native contacts, e.g. CA,CB,C,N,O  (%default)")
	op.add_option('-n','--chainSep', action="store", type="int", dest="chainSep",help="min number of residues between contacting atoms  (%default)")
	op.add_option('-d','--distCutoff', action="store", type="float", dest="distCutoff",help="distance cutoff for defining contacts (%default Ang)")
	op.add_option('-m','--wellDepth', action="store", type="float", dest="wellDepth",help="well depth at native distance (%default)")
	op.add_option('-w','--wellWidth', action="store", type="float", dest="wellWidth",help="GAUSS/DUALGAUSS: width of well (%default)")
	op.add_option('-r','--repulDist', action="store", type="float", dest="repulDist",help="GAUSS/DUALGAUSS: excluded volume radius.  0=no repulsion; -x=radius is min distance per contact - x (FMULTIGAUSS)  (%default)")
	op.add_option('-p','--repulSteep', action="store", type="float", dest="repulSteep",help="GAUSS/DUALGAUSS: steepness of repulsion (%default)")
	op.add_option('-e','--ensembleCutoff', action="store", type="float", dest="ensembleCutoff",help="when PDB file contains multiple models, retain contacts present in this fraction of models (%default)")
	op.add_option('-N','--normEnergy', action="store",type="float", dest="normEnergy",help="normalize contact energies to given amount per structure. Overrides options m and M. (%default=don't normalize)")
	op.add_option('--gmxpath', action="store", type="string", dest="gmxpath",help="path to Gromacs executables (SHADOW)")
	op.add_option('--scmpath', action="store", type="string", dest="scmpath",help="path to SCM.jar  (SHADOW)")
	op.add_option('-L','--label',action='store', type="string",dest="label",help="string that will be appended to the keyword 'DistanceRestraint' in the Profasi energy output")
	op.add_option('-S','--seqFromModel',action='store', type="int",dest="seqFrom",help="model number to select AA sequence for contacts. By default: take first model")
	op.add_option('-x','--nonNative',action='store_true',dest="includeNonNative",help="include non-native contacts")
	op.add_option('-R','--range',action='store',type="string",dest="selectedRange",help="specify residues ranges as strings,e.g. '1-20,29,30,42-56' add 'x' at the beginning to exclusively only consider these residues, otherwise contacts involving only one of these residues will be considered")
	op.add_option('-M','--smoothGaussian',action='store_true',dest="smoothGaussian",help='set function to epsilon if between min and max distance')
	op.set_defaults(inPDBFilename="",
					contactType="CUTOFF",
					contactAtoms="CA",
					chainSep=3,
					distCutoff=6.0,
					wellDepth=1.0,
					wellWidth=0.5,
					repulDist=0.0,
					repulSteep=1.0,
					ensembleCutoff=0.8,
					normEnergy=0,
					gmxpath="",
					scmpath="",
					label="",
					seqFrom=0,
					includeNonNative=False,
					selectedRange="",
					smoothGaussian=False)
	
	opt, args = op.parse_args()
								
	
	inPDBFilename = opt.inPDBFilename
	contactType = opt.contactType
	contactAtoms = opt.contactAtoms
	chainSep = opt.chainSep
	distCutoff = opt.distCutoff
	wellDepth = opt.wellDepth
	wellWidth = opt.wellWidth
	repulDist = opt.repulDist
	repulSteep = opt.repulSteep
	ensembleCutoff = opt.ensembleCutoff
	normEnergy = opt.normEnergy
	gmxpath = opt.gmxpath
	scmpath = opt.scmpath
	label = opt.label
	seqFrom = opt.seqFrom
	includeNonNative = opt.includeNonNative
	selectedRange = opt.selectedRange.strip("'").strip('"')
	smoothGaussian = opt.smoothGaussian
	
	assert contactType in ["CUTOFF","SHADOW","HYBRID","RESIDUE"]
	assert distCutoff > 1.0, "cutoff in Angstrom!"
	
	restraintType = "native_restraints"
	assert restraintType in ["distance_restraints", "native_restraints"]
	
	contactPotential = "FMULTIGAUSS"
	if smoothGaussian:
		contactPotential = "FMULTIGSMOOTH"
	assert contactPotential in ["FMULTIGAUSS","FMULTIGSMOOTH"]
	
	selection = []
	selectExclusive = False
	if selectedRange!="":
		
		tmp = selectedRange.strip().split(",")
		if tmp[0][0]=="x":
			tmp[0] = tmp[0].strip("x")
			selectExclusive = True
		for t1 in tmp:
			if not "-" in t1:
				selection.append(int(t1))
			else:
				tmp2 = t1.strip().split("-")
				assert len(tmp2)==2
				start = int(tmp2[0])
				end = int(tmp2[1])
				for i in range(start,end+1):
					selection.append(i)
		print "Your selection of amino acid positions: ", selectedRange, selection
	
		
	consStatement = "_cons%s"%str(ensembleCutoff)
	
	
	removeShared = None
	
	contactAtoms = contactAtoms.split(",")
	
	if contactType == "SHADOW" and contactAtoms != ['CA']:
		print "SHADOW only allows CA atoms!"
		contactAtoms = ['CA']
	
	if label=="":
		label = "_"+inPDBFilename[:-4]
	
	print "Atoms to consider for CUTOFF contacts:",contactAtoms
	
	assert all([ a in BB_atoms for a in contactAtoms])
	
	if contactType in ["HYBRID","SHADOW"]:
		print "Also considering CA atoms for SHADOW contacts."
	
	if normEnergy > 0: print "Normalizing contact energy per native structure to:", normEnergy
	
	assert inPDBFilename != ""
	
	
	atoms_models = []
	models_contact_lists = []
	models_scenergy_lists = []
	
	#atoms_models = extractModelAtoms(inPDBFilename)
	models = splitModels(inPDBFilename)
	
	print "\n\n%i model(s) found in %s."%(len(models), inPDBFilename)
	
	for m in models:
		contacts = []
		con_c=[]
		con_s=[]
		
		atoms = atomsFromPDB(m)
		atoms_models.append(atoms)
		
		print "num atoms",len(atoms), m
		
		if contactType=="CUTOFF":
			m_contacts = getCutoffContacts(atoms, atomtypes=contactAtoms, cutoff=distCutoff, chainDist=chainSep, verbose=False, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		elif contactType=="SHADOW":
			m_contacts = getShadowContacts(pdbfilename=m, outfilename=inPDBFilename[:-4]+".shadow", cutoff=distCutoff, chainDist=3, gmxpath=gmxpath, scmpath=scmpath, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		elif contactType=="HYBRID":
			con_c = getCutoffContacts(atoms, atomtypes=contactAtoms, cutoff=distCutoff, chainDist=chainSep, verbose=False, res2SS = None, selection=selection, selectExclusive=selectExclusive)
			con_s = getShadowContacts(pdbfilename=m, outfilename=inPDBFilename[:-4]+".shadow", cutoff=distCutoff, chainDist=3, gmxpath=gmxpath, scmpath=scmpath, res2SS = None, selection=selection, selectExclusive=selectExclusive)
			m_contacts = mixCutoffShadow(con_c, con_s)
		elif contactType=="RESIDUE":
			m_contacts = getResidueContacts(atoms, cutoff=distCutoff, chainDist=chainSep,verbose=False, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		else:
			print "unknown contact type: ",contactType
			sys.exit(1)
		
		
		models_contact_lists.append(m_contacts)
		
		
		
		
	
	
	
	
	print "%s ensemble contacts"%contactPotential
	finalcontacts = filterSuperEnsembleContacts(models_contact_lists, atoms_models, contactAtoms, ensembleCutoff, contactType)
	#print finalcontacts
	
	
	for m in models:
		energies_m =[]# getSCenergies4Contacts(finalcontacts, m, atoms_models[seqFrom])
		#print "Total SC energy: ", sum([e[-1] for e in energies_m])
		
		models_scenergy_lists.append( energies_m  )
	
	pbond_legend = None
	
	
	outatoms = "".join([i for i in contactAtoms])
	
	
	
	width_label = "_w"+str(round(wellWidth,1))
	repul_label = "_R"+str(round(repulDist,1))
	nonnat_label = "_x"+str(includeNonNative)
	
	
	outfilename = "%s_%s_%s_c%s_n%s%s%s_%s%s%s_%s"%(inPDBFilename.split("/")[-1],contactPotential,contactType,str(distCutoff),str(chainSep),width_label,repul_label,outatoms,consStatement,nonnat_label,selectedRange)
	
	if contactType=="HYBRID" and not 'CA' in contactAtoms:
		contactAtoms.append('CA')
	
	
	
	
	print "Sequence from model",seqFrom
	
	outfilename = writeProfasiContactsFMULTIGAUSS(finalcontacts, atoms_models[seqFrom], contactAtoms, outfilename, repulDist, repulSteep, width=wellWidth, depth=wellDepth, norm=normEnergy,label=label,includeNonNative=includeNonNative,chainSep=chainSep,smoothGaussian=smoothGaussian,restraintType=restraintType)
	
	writeChimeraPseudobonds(finalcontacts, outfilename[:-4]+".pseudobonds",legend=pbond_legend)
	
	print "wrote Profasi XML contacts:",outfilename
	print "wrote Chimera pseudobonds:",outfilename[:-4]+".pseudobonds"
	#if not potentialType=="FMULTIGAUSS":
	#	verifyProfasiXML(outfilename, bestModels, contactAtoms)
	if True:
		plotMultiGaussianContacts(finalcontacts,1.0,0.5,outfilename[:-4]+"_multiwell.png",smoothGaussian=smoothGaussian)
	if True:
		nres = getLength(atoms)
		plotContactMap(nres, finalcontacts, inPDBFilename, distCutoff, chainSep, atoms, contactType, fname=outfilename[:-4]+".png")
	
	
	
	
	