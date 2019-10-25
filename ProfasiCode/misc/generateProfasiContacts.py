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
			if os.path.exists(modelfilename):
				model_exists = True
			
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
	
	

def filterSuperEnsembleContacts(models_contact_lists, atoms, atypes, ensembleCutoff, contactType, verbose=False):
	print "super ensemble contacts"
	contactsDict = {}
	
	atomtypes = copy.copy(atypes)
	
	if contactType == "HYBRID" and not "CA" in atomtypes:
		atomtypes  += ['CA']
	
	nModels = len(models_contact_lists)
	contactsPerModel = [[] for i in xrange(nModels)]
	
	ref_atoms = atoms[0]
	
	for m in xrange(nModels):
		model_contacts = models_contact_lists[m]
		#print "%i contacts in model %s"%(len(model_contacts),m)
		for c in model_contacts:
			a1 = c[0]
			a2 = c[1]
			d  = c[2]
			modelatoms = atoms[m]
			assert modelatoms[a1]["aname"] in atypes, modelatoms[a1]
			assert modelatoms[a2]["aname"] in atypes, modelatoms[a2]
			#contactsPerModel[m].append( [a1,a2,d,modelatoms] )
			
			
			ref_a1,ref_a2 = translateAtoms(a1,a2,atoms[m],ref_atoms)
			
			if ref_a1==None:
				print "Could not find atom %i from model %i in model 0! Skipping this contact!"%(a1,m)
				continue
			if ref_a2==None:
				print "Could not find atom %i from model %i in model 0! Skipping this contact!"%(a2,m)
				continue
			
			if (ref_a1,ref_a2) in contactsDict:
				contactsDict[(ref_a1,ref_a2)].append(d)
			else:
				contactsDict[(ref_a1,ref_a2)] = [d]
	
	# which contacts are found in all models? what are the means and standard deviations of distances?

	counter = 0
	avgStdDistsShared = []
	
	finalContacts = []
	dists = []
	means = []
	stddevs = []
	
	global_min_dist = 10000
	global_max_dist = -1
	min_range = 100000
	max_range = -1
	
	for i in sorted(contactsDict.keys()):
	
		contactCount = len(contactsDict[i]) 
		assert contactCount <= nModels
		
		if contactCount / float(nModels) >= ensembleCutoff: # in how many models was this contact found ? Is the fraction above threshold?
			counter += 1
			
			avgDist = numpy.mean(contactsDict[i])
			stdDist = numpy.std(contactsDict[i])
			
			
			min_d = min(contactsDict[i])
			max_d = max(contactsDict[i])
			range_d = max_d - min_d
			
			if min_d < global_min_dist: global_min_dist = min_d
			if max_d > global_max_dist: global_max_dist = max_d
			if range_d < min_range: min_range = range_d
			if range_d > max_range: max_range = range_d
			
			print "atom pair:", i, " -- ", str((ref_atoms[i[0]]["aname"], (ref_atoms[i[1]]["aname"]))),"  ", contactCount, "models with this contact. Avg,Std:", avgDist, stdDist, "Min,Max:",min_d,max_d," Range:",range_d#, 100*stdDist/avgDist, "%"
			
			
			
			if ref_atoms[i[0]]["aname"] in atomtypes and ref_atoms[i[1]]["aname"] in atomtypes:
				#assert atoms[i[0]]["aname"] in atomtypes
				#assert atoms[i[1]]["aname"] in atomtypes
				#print contactsDict[i]
			
				finalContacts.append([i[0], i[1], avgDist])
				
				dists.append(contactsDict[i])
				means.append(avgDist)
				stddevs.append(stdDist)
			
				avgStdDistsShared.append([i[0], i[1], avgDist, stdDist])
			else:
				print "atom is not allowed for contact!",ref_atoms[i[0]]["aname"],ref_atoms[i[1]]["aname"]
	
	print "Min and Max distances within all contacts:", global_min_dist, global_max_dist
	print "Min and Max distance ranges:", min_range,max_range
	
	print "%i contacts above threshold"%(counter)
	
	
	
	return finalContacts, ref_atoms, dists, means, stddevs

def filterEnsembleContacts(models_contact_lists, atoms, atypes, ensembleCutoff, contactType, verbose=False):
	contactsDict = {}
	
	atomtypes = copy.copy(atypes)
	
	if contactType == "HYBRID" and not "CA" in atomtypes:
		atomtypes  += ['CA']
	
	nModels = len(models_contact_lists)
	contactsPerModel = [[] for i in xrange(nModels)]
	
	for m in xrange(nModels):
		model_contacts = models_contact_lists[m]
		#print "%i contacts in model %s"%(len(model_contacts),m)
		for c in model_contacts:
			a1 = c[0]
			a2 = c[1]
			d  = c[2]
			
			contactsPerModel[m].append( [a1,a2,d] )
			
			if (a1,a2) in contactsDict:
				contactsDict[(a1,a2)].append(d)
			else:
				contactsDict[(a1,a2)] = [d]
	
	# which contacts are found in all models? what are the means and standard deviations of distances?

	counter = 0
	avgStdDistsShared = []
	
	for i in sorted(contactsDict.keys()):
	
		contactCount = len(contactsDict[i]) 
		assert contactCount <= nModels
		
		if contactCount / float(nModels) >= ensembleCutoff: # in how many models was this contact found ? Is the fraction above threshold?
			counter += 1
			avgDist = numpy.mean(contactsDict[i])
			stdDist = numpy.std(contactsDict[i])
			if verbose: print "atom pair:", i, " -- ", str((atoms[i[0]]["aname"], (atoms[i[1]]["aname"]))),"  ", contactCount, "models with this contact. Avg,Std:", avgDist, stdDist#, 100*stdDist/avgDist, "%"
			assert atoms[i[0]]["aname"] in atomtypes
			assert atoms[i[1]]["aname"] in atomtypes
			#print contactsDict[i]
			avgStdDistsShared.append([i[0], i[1], avgDist, stdDist])
	
	print "%i contacts above threshold"%(counter)
	
	minDifference = 100000
	minDiffModel = None
	for i in xrange(len(contactsPerModel)):
		contacts = contactsPerModel[i] # all contacts from ONE model
		summedSQdiff = 0
		for j in xrange(len(avgStdDistsShared)):
			shared_contact = avgStdDistsShared[j][:2]
			avgDist = avgStdDistsShared[j][2]
			stdDist = avgStdDistsShared[j][3]
		
			for k in contacts:
				#print k[:2] , shared_contact
				if k[:2] == shared_contact:
					diff = k[2] - avgDist
					sqdiff = math.pow(diff,2)
					summedSQdiff += sqdiff
					#print k[:2], diff
		if verbose: print "model ",i+1," -- summed sq diff:", summedSQdiff
		if summedSQdiff < minDifference:
			minDifference = summedSQdiff
			minDiffModel = i

	if verbose: print "min",minDifference,minDiffModel + 1

	print "Selecting Model with smallest squared difference of contact distances compared to ensemble average! --> Model ",minDiffModel+1

	finalContacts = []
	
	
	# which contacts from the most representative model are shared by the required fraction of models ?
	contacts = contactsPerModel[minDiffModel]
	print len(contacts), "contacts originally in model",minDiffModel+1
	for c in contacts:
		#print c
		for j in xrange(len(avgStdDistsShared)):
			shared_contact = avgStdDistsShared[j][:2]
			if c[:2] == shared_contact:
				finalContacts.append(c)
			
	
	print "No. of contacts in model %i found in %f of %i models: %i"%(minDiffModel + 1,ensembleCutoff,nModels,len(finalContacts))
	
	return finalContacts, minDiffModel
	

def prf_convert(infile, outfile, prf_loc):
	exec_string = "%sprf_convert %s -o %s"%(prf_loc, infile, outfile)
	sub = subprocess.Popen(exec_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out = sub.stdout.read()
	err = sub.stderr.read()
	print err



def regularize(inpdb, prf_loc):
	exec_string = "%sregularize %s"%(prf_loc, inpdb)
	sub = subprocess.Popen(exec_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out = sub.stdout.read()
	err = sub.stderr.read()
	print err
	
	prf_convert("min_rmsd.xml", inpdb[:-4]+"_min_rmsd.pdb", prf_loc)
	prf_convert("min_rmsd.xml", inpdb[:-4]+"_min_rmsd.tconf", prf_loc)
	prf_convert("min_etot.xml", inpdb[:-4]+"_min_etot.pdb", prf_loc)
	prf_convert("min_etot.xml", inpdb[:-4]+"_min_etot.tconf", prf_loc)



def getCutoffContacts(pdbfilename, atomtypes, cutoff, chainDist, verbose=False, res2SS = None):
	atoms = atomsFromPDB(pdbfilename)
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
		if res2SS != None and res2SS[ri][0] == 'X':
			#print "skip ri",ri,res2SS[ri]
			continue

		for j in xrange( len( contactAtoms ) ):
			aj = contactAtoms[j]
			rj = atoms[ aj ][ "rnum" ]
			if res2SS != None and res2SS[rj][0] == 'X':
				#print "skip rj",rj,res2SS[rj]
				continue
			
			if ai < aj:
				if abs(rj-ri) >= chainDist:
					d = atomDist( ai, aj, atoms )
					#print d
					if d <= cutoff:
						contacts.append( [ai, aj, d] )
						if verbose: print [ai, aj, d], ri,rj
	#print contacts
	
	print "Found %i CUTOFF contacts in %s between %s atoms at chain separation=%i and cutoff distance = %f"%(len(contacts), pdbfilename, str(atomtypes),chainDist, cutoff)
	
	return contacts

def getResidueContacts(pdbfilename, cutoff, chainDist, atoms, verbose=False, res2SS = None):
	from Bio.PDB import *
	
	structure = PDBParser.get_structure(PDBParser(),pdbfilename,pdbfilename)

	cutoff = cutoff

	chainsep = chainDist

	model = structure[0]

	reslist = Selection.unfold_entities(structure, 'R')

	contacts = []
	
	
	for i in xrange(len(reslist)):
		ri = reslist[i]
		for j in xrange(len(reslist)):
			rj = reslist[j]
		
			if i<j-chainsep:
				mindist = 10000
				mdpair = None
				cadist = None
				for ai in ri:
					if 'H' in ai.get_id(): 
						#print ai.get_id()
						continue
					for aj in rj:
						if 'H' in aj.get_id(): continue
						d = abs(aj-ai)
					
						if cadist == None:
							if ai.get_id() == 'CA' and ai.get_id() == 'CA':
								cadist = d
						if d < mindist:
							mindist=d
							mdpair = (ai,aj)
				if mindist <= cutoff:
					a1 = getAtomIDforResidueContact(i+1,'CA',atoms)
					a2 = getAtomIDforResidueContact(j+1,'CA',atoms)
					
					contacts.append( [a1,a2,cadist] )
					print [ri,rj,mindist,cadist], mdpair[0].get_id(), mdpair[1].get_id()	,a1,a2
	
	
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
	


def getShadowContacts(pdbfilename, outfilename, cutoff=6.0, chainDist=3, gmxpath="", scmpath="", res2SS = None, verbose = False):
	
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

def mixContactsBistable(contacts1, contacts2, removeShared=True, legend=None):

	if not removeShared:
		if legend != None:
			legend.extend(["s1" for i in xrange(len(contacts1))])
			legend.extend(["s2" for i in xrange(len(contacts2))])
		return contacts1 + contacts2
	else:
		allContacts = []
		pairs1 = {}
		for c in contacts1: 
			p = tuple(c[:2])
			pairs1[p]=True

		pairs2 = {}
		for c in contacts2:
			p = tuple(c[:2])
			pairs2[p]=True

		for c in contacts1: 
			p = tuple(c[:2])
			if not p in pairs2:
				legend.append("s1")
				allContacts.append(c)
		for c in contacts2: 
			p = tuple(c[:2])
			if not p in pairs1:
				legend.append("s2")
				allContacts.append(c)

		return allContacts

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

def writeProfasiContactsLJ(contacts, atoms, outfilename, depth, norm=0,label=""):
	out = open(outfilename,"w")
	out.write("""
<distance_restraints rid="%s">  
<formatted_data>"
    <format name=\"restraint\" type=\"$3\">
      <atom1>$1</atom1>"
      <atom2>$2</atom2>
      <parameters>
        <epsilon>$4</epsilon>
        <minimum>$5</minimum>
      </parameters>
    </format>
<data>
"""%(label))
	n = len(contacts)
	if norm > 0:
		depth = norm/float(n)
	
	for c in contacts:
		a1 = c[0]
		a2 = c[1]
		d = c[2]
		
		atomstring1 = profasiAtomString(a1, atoms)
		atomstring2 = profasiAtomString(a2, atoms)
		out.write( "  %s  %s  LJ  %s  %s\n"%(atomstring1, atomstring2, str(depth), str(d)) )
	out.write("""
    </data>
  </formatted_data>
</distance_restraints>
""")


def writeProfasiContactsGAUSS(contacts, atoms, outfilename, radius, steepness, width, depth, norm=0,label=""):
	out = open(outfilename,"w")
	out.write("""
<distance_restraints rid="%s">  
<formatted_data>"
    <format name=\"restraint\" type=\"$3\">
      <atom1>$1</atom1>"
      <atom2>$2</atom2>
      <parameters>
        <minimum>$4</minimum>
        <radius>$5</radius>
        <steepness>$6</steepness>
        <width>$7</width>
        <depth>$8</depth>
      </parameters>
    </format>
<data>
"""%(label))
	n = len(contacts)
	if norm > 0:
		depth = norm/float(n)
		print "New well depth:",depth
	
	for c in contacts:
		a1 = c[0]
		a2 = c[1]
		d = c[2]
		
		
		
		atomstring1 = profasiAtomString(a1, atoms)
		atomstring2 = profasiAtomString(a2, atoms)
		out.write( "  %s  %s  GAUSS  %s  %s  %s %s %s\n"%(atomstring1, atomstring2, str(d), str(radius), str(steepness), str(width), str(depth) ) )
	
	out.write("""
    </data>
  </formatted_data>
</distance_restraints>
""")

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

def writeProfasiContactsDUALGAUSS(contacts1, contacts2, atoms1, atoms2, atomtypes, outfilename, radius, steepness, width, depth, width2, depth2, norm=0,verbose=False,label=""):
	out = open(outfilename,"w")
	out.write("""
<distance_restraints rid="%s">  
<formatted_data>"
    <format name=\"restraint\" type=\"$3\">
      <atom1>$1</atom1>"
      <atom2>$2</atom2>
      <parameters>
        <minimum>$4</minimum>
        <minimum2>$5</minimum2>
        <radius>$6</radius>
        <steepness>$7</steepness>
        <width>$8</width>
        <depth>$9</depth>
        <width2>$10</width2>
        <depth2>$11</depth2>
      </parameters>
    </format>
<data>
"""%(label))
	
	#a2a = atoms2atoms(atoms2,atoms1,atomtypes) # in case atoms in the two pdb structures were numbered differently !!!
	
	potential = "DUALGAUSS"
	parameters = [None, None, radius, steepness, width, depth, width2, depth2]
	
	energies1 = [] # energies when structure is structure 1
	energies2 = [] # or 2
	energies1_shared = []
	energies2_shared = []
	
	allcontacts = []
	shared = []
	unique1 = []
	unique2 = []
	
	for c1 in contacts1:
		c1_a1 = c1[0]
		c1_a2 = c1[1]
		d1 = c1[2]
		
		
		isShared = False
		for c2 in contacts2:
			c2_a1 = c2[0]
			c2_a2 = c2[1]
			d2 = c2[2]
			
			if (c1_a1, c1_a2) == (c2_a1, c2_a2):
				isShared = True
				if verbose: print "shared contacts:",c1,c2
				shared.append((c1,c2))
				
				
				#allcontacts.append((c1,c2))
				#if legend != None: legend.append("s12")
		
		if not isShared:
			#allcontacts.append(c1)
			#if legend != None: legend.append("s1")
			unique1.append(c1)
			
			
				
	for c2 in contacts2:
		c2_a1 = c2[0]
		c2_a2 = c2[1]
		d2 = c2[2]
		
		found_in_shared = False
		for s in shared:
			s_a1 = s[0][0]
			s_a2 = s[0][1]
			if (s_a1, s_a2) == (c2_a1, c2_a2):
				found_in_shared = True
				break
		if not found_in_shared:
			#allcontacts.append(c2)
			#if legend != None: legend.append("s2")
			unique2.append(c2)
			
	
	n1 = len(unique1)+len(shared)
	n2 = len(unique2)+len(shared)
	
	print "Contacts in native state of 1 or 2:",n1,n2
	
	print "Old well depths:",depth,depth2
	if norm > 0:
		idealdepth = norm/float(n1)
		idealdepth2 = norm/float(n2)
		
		print "New well depths:",idealdepth,idealdepth2
	else:
		idealdepth = depth
		idealdepth2 = depth2
	
	for ui in xrange(len(unique1)):
		u1 = unique1[ui]
		a1 = u1[0]
		a2 = u1[1]
		d1 = u1[2]
		e1 = getEnergy(potential, [d1, d1, radius, steepness, width, idealdepth, 0, 0], d1 )
		
		energies1.append(e1)
	for ui in xrange(len(unique2)):
		u2 = unique2[ui]
		a1 = u2[0]
		a2 = u2[1]
		d2 = u2[2]
		e2 = getEnergy(potential, [d2, d2, radius, steepness, 0, 0, width2, idealdepth2], d2 )
		
		energies2.append(e2)
	for si in xrange(len(shared)):
		s = shared[si]
		assert s[0][0] == s[1][0]
		assert s[0][1] == s[1][1]
		a1 = s[0][0]
		a2 = s[0][1]
		
		d1 = s[0][2]
		d2 = s[1][2]
		
		e1 = getEnergy(potential, [d1, d2, radius, steepness, width, idealdepth, width2, idealdepth2], d1 )
		e2 = getEnergy(potential, [d1, d2, radius, steepness, width, idealdepth, width2, idealdepth2], d2 )
		
		energies1_shared.append(e1)
		energies2_shared.append(e2)
	
	
	
	
	
	print sum(energies1)
	print sum(energies2)
	print sum(energies1_shared)
	print sum(energies2_shared)
	print sum(energies1) + sum(energies1_shared)
	print sum(energies2) + sum(energies2_shared)
	
	corr1 = []
	corr2 = []
	
	#contactsBoth.extend(allcontacts)
	
	for ui in xrange(len(unique1)):
		u1 = unique1[ui]
		
		a1 = u1[0]
		a2 = u1[1]
		d = u1[2]
		
		e1 = energies1[ui]
		
		correction1 = e1 - idealdepth
		
		corr_depth = abs(depth + correction1)
		corr1.append(corr_depth)
		
		atomstring1 = profasiAtomString(a1, atoms1)
		atomstring2 = profasiAtomString(a2, atoms1)
		out.write( "  %s  %s  DUALGAUSS  %s  %s  %s  %s  %s  %s  %s  %s\n"%(atomstring1, atomstring2, str(d), str(d), str(radius), str(steepness), str(width), str(idealdepth), str(0), str(0)) )
	
	for ui in xrange(len(unique2)):
		u2 = unique2[ui]
		a1 = u2[0]
		a2 = u2[1]
		d = u2[2]
		
		e2 = energies2[ui]
		
		correction2 = e2 - idealdepth2
		
		corr_depth2 = abs(depth2 + correction2)
		
		corr2.append(corr_depth2)
		
		atomstring1 = profasiAtomString(a1, atoms2)
		atomstring2 = profasiAtomString(a2, atoms2)
		out.write( "  %s  %s  DUALGAUSS  %s  %s  %s  %s  %s  %s  %s  %s\n"%(atomstring1, atomstring2, str(d), str(d), str(radius), str(steepness), str(0), str(0), str(width2), str(idealdepth2)) )
	
	for si in xrange(len(shared)):
		s = shared[si]
		assert s[0][0] == s[1][0]
		assert s[0][1] == s[1][1]
		a1 = s[0][0]
		a2 = s[0][1]
		
		d1 = s[0][2]
		d2 = s[1][2]
		
		e1 = energies1_shared[si]
		e2 = energies2_shared[si]
		
		correction1 = e1 + idealdepth
		correction2 = e2 + idealdepth2
		print correction1, correction2
		
		corr_depth = idealdepth# - correction1
		corr_depth2 = idealdepth2# - correction2
		print corr_depth,corr_depth2
		
		corr1.append(corr_depth)
		corr2.append(corr_depth2)
		
		e1_corr = getEnergy(potential, [d1, d2, radius, steepness, width, corr_depth, width2, corr_depth2], d1 )
		e2_corr = getEnergy(potential, [d1, d2, radius, steepness, width, corr_depth, width2, corr_depth2], d2 )
		print e1, e1_corr, e2, e2_corr
		
		atomstring1 = profasiAtomString(a1, atoms1)
		atomstring2 = profasiAtomString(a2, atoms1)
		out.write( "  %s  %s  DUALGAUSS  %s  %s  %s  %s  %s  %s  %s  %s\n"%(atomstring1, atomstring2, str(d1), str(d2), str(radius), str(steepness), str(width), str(corr_depth), str(width2), str(corr_depth2)) ) 

	out.write("""
    </data>
  </formatted_data>
</distance_restraints>
""")

	print sum(corr1)
	print sum(corr2)
	
	print "DUALGAUSS contacts - unique to 1: %d, unique to 2: %d, shared by 1 and 2: %d (total: %d  +18=%d)"%(len(unique1), len(unique2), len(shared), len(unique1)+len(unique2)+len(shared), len(unique1)+len(unique2)+len(shared)+len(shared))
	
	print "ENERGY = U1 + U2 = %f + %f = %f"%(len(unique1)*depth+len(shared)*depth, len(unique2)*depth2+len(shared)*depth2, len(unique1)*depth+len(shared)*depth + len(unique2)*depth2+len(shared)*depth2)


def writeProfasiContactsFDUALGAUSS(contacts1, contacts2, atoms1, atoms2, atomtypes, outfilename, radius, steepness, width, depth, width2,verbose=False,label=""):
	out = open(outfilename,"w")
	out.write("""
<distance_restraints rid="%s">  
<formatted_data>"
    <format name=\"restraint\" type=\"$3\">
      <atom1>$1</atom1>"
      <atom2>$2</atom2>
      <parameters>
        <minimum>$4</minimum>
        <minimum2>$5</minimum2>
        <radius>$6</radius>
        <steepness>$7</steepness>
        <width>$8</width>
        <depth>$9</depth>
        <width2>$10</width2>
        
      </parameters>
    </format>
<data>
"""%(label))
	
	#a2a = atoms2atoms(atoms2,atoms1,atomtypes) # in case atoms in the two pdb structures were numbered differently !!!
	
	
	
	
	allcontacts = []
	shared = []
	unique1 = []
	unique2 = []
	
	for c1 in contacts1:
		c1_a1 = c1[0]
		c1_a2 = c1[1]
		isShared = False
		for c2 in contacts2:
			c2_a1 = c2[0]
			c2_a2 = c2[1]
			
			if (c1_a1, c1_a2) == (c2_a1, c2_a2):
				isShared = True
				if verbose: print "shared contacts:",c1,c2
				shared.append((c1,c2))
				#allcontacts.append((c1,c2))
				#if legend != None: legend.append("s12")
		
		if not isShared:
			#allcontacts.append(c1)
			#if legend != None: legend.append("s1")
			unique1.append(c1)
				
	for c2 in contacts2:
		c2_a1 = c2[0]
		c2_a2 = c2[1]
		
		found_in_shared = False
		for s in shared:
			s_a1 = s[0][0]
			s_a2 = s[0][1]
			if (s_a1, s_a2) == (c2_a1, c2_a2):
				found_in_shared = True
				break
		if not found_in_shared:
			#allcontacts.append(c2)
			#if legend != None: legend.append("s2")
			unique2.append(c2)
	
	
	
	n1 = len(unique1)+ len(unique2)+len(shared)
	
	
	print "Contacts in native state of 1 or 2:",n1
	
	
	
	
	#contactsBoth.extend(allcontacts)
	
	for u in unique1:
		a1 = u[0]
		a2 = u[1]
		d = u[2]
		
		
		
		atomstring1 = profasiAtomString(a1, atoms1)
		atomstring2 = profasiAtomString(a2, atoms1)
		out.write( "  %s  %s  FDUALGAUSS  %s  %s  %s  %s  %s  %s  %s\n"%(atomstring1, atomstring2, str(d), str(d), str(radius), str(steepness), str(width), str(depth), str(0)) )
	
	for u in unique2:
		a1 = u[0]
		a2 = u[1]
		d = u[2]
		
		
		
		atomstring1 = profasiAtomString(a1, atoms2)
		atomstring2 = profasiAtomString(a2, atoms2)
		out.write( "  %s  %s  FDUALGAUSS  %s  %s  %s  %s  %s  %s  %s\n"%(atomstring1, atomstring2, str(d), str(d), str(radius), str(steepness), str(0), str(depth), str(width2)) )
	
	for s in shared:
		assert s[0][0] == s[1][0]
		assert s[0][1] == s[1][1]
		a1 = s[0][0]
		a2 = s[0][1]
		
		d1 = s[0][2]
		d2 = s[1][2]
		
		atomstring1 = profasiAtomString(a1, atoms1)
		atomstring2 = profasiAtomString(a2, atoms1)
		out.write( "  %s  %s  FDUALGAUSS  %s  %s  %s  %s  %s  %s  %s\n"%(atomstring1, atomstring2, str(d1), str(d2), str(radius), str(steepness), str(width), str(depth), str(width2)) ) 

	out.write("""
    </data>
  </formatted_data>
</distance_restraints>
""")
	print "FDUALGAUSS contacts - unique to 1: %d, unique to 2: %d, shared by 1 and 2: %d (total: %d  +18=%d)"%(len(unique1), len(unique2), len(shared), len(unique1)+len(unique2)+len(shared), len(unique1)+len(unique2)+len(shared)+len(shared))
	
	print "ENERGY = U1 + U2 = %f + %f = %f"%(len(unique1)*depth+len(shared)*depth, len(unique2)*depth+len(shared)*depth, len(unique1)*depth+len(shared)*depth + len(unique2)*depth+len(shared)*depth)

def writeProfasiContactsFMULTIGAUSS(contacts, dists, atoms, outfilename, radius, steepness, width, depth, norm=0,label=""):
	out = open(outfilename,"w")
	out.write("""
<distance_restraints rid="%s">  
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
<data>
"""%(label))
	n = len(contacts)
	if norm > 0:
		depth = norm/float(n)
		print "New well depth:",depth

	for ci in xrange(len(contacts)):
		c = contacts[ci]
		a1 = c[0]
		a2 = c[1]
		d = c[2]
		dist = dists[ci]
		
		mindist = min(dist)
		if radius < 0:
			final_radius = mindist + radius
		else:
			final_radius = radius
		
		atomstring1 = profasiAtomString(a1, atoms)
		atomstring2 = profasiAtomString(a2, atoms)
		out.write( "  %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n"%(atomstring1, atomstring2, ",".join([str(i) for i in dist]), str(final_radius), str(steepness), str(width), str(depth) ) )

	out.write("""
    </data>
  </formatted_data>
</distance_restraints>
""")
	print "Expected energy:",n*depth
	


def writeChimeraPseudobonds(contacts, atoms, outfilename, legend=None):

	outfile = open(outfilename, 'w')
	
	if legend != None: assert len(contacts) == len(legend),str(len(contacts))+" "+str(len(legend))
	
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		a1 = c[0]
		a2 = c[1]
		
		if len(c)>2:
			d  = c[2]
			label = str(round(d,2))
		else:
			label = ""
		
		if potentialType == "FMULTIGAUSS":
			label = str(round(min(dists[ci]),2)) + "-" + str(round(max(dists[ci]),2))
		
		if legend==None:
			color = "blue"
			
		else:
			lstring = legend[ci]
			if lstring == "s1":
				color = "blue"
			elif lstring == "s2":
				color = "red"
			else:
				color = "magenta"
		pseudobond = "#0:%s@%s #0:%s@%s %s %s"%(str(atoms[a1]["rnum"]), atoms[a1]["aname"].lower(), str(atoms[a2]["rnum"]), atoms[a2]["aname"].lower(), color, label)
		outfile.write(pseudobond + "\n")

	outfile.close()




def verifyProfasiXML(outfilename, pdb_files, atomtypes, tolerance=0.01):
	
	print "verifying contacts..."
	xml_lines = []
	start=False
	distances = []
	distances2 = []
	multidistances = []
	
	xml_fulllines = []
	for line in open(outfilename).readlines():
		if "<data>" in line:
			#print "found data"
			start=True
		elif "</data>" in line:
			#print "finished reading data"
			start = False
		elif start:
			tmp = line.strip().split()
			if len(tmp) >=3: 
				xml_fulllines.append(tmp)
				xml_lines.append(tmp[:2])
				potential = tmp[2]
				if potential=="FMULTIGAUSS":
					
					
					multidistances.append([float(t) for t in tmp[3].split(",")])
				else:
					distances.append(float(tmp[3]))
				if potential=="DUALGAUSS" or potential=="FDUALGAUSS":
					#print line
					distances2.append(float(tmp[4]))
	
	print "Found %i contacts in %s"%(len(xml_lines),outfilename)
	
	assert len(xml_lines)==len(distances)
	
	distances_list = [distances]
	if len(pdb_files)==2:
		assert len(distances2)==len(xml_lines)
		distances_list.append(distances2)
	
	
	
	avgAbsDiff = []
	countDiff = 0
	avgAbsDiff2 = []
	countDiff2 = 0
	nFound=0
	nFound2=0
	
	totalEnergy1 = []
	totalEnergy2 = []
	
	
	
	
	
	for x in xrange(len(xml_lines)):
		i = xml_lines[x]
		t = xml_fulllines[x][2]
		parms = [float(xx) if potential!="FMULTIGAUSS" else xx for xx in xml_fulllines[x][3:] ]
		print xml_fulllines[x]
		
		tmp1=i[0].split("/")
		tmp2=i[1].split("/")
		
		rnum1 = int(tmp1[1])+1 # PROFASI COUNTS RESIDUE NUMBERS FROM 0 !!!
		rname1 = tmp1[2]
		aname1 = tmp1[3].strip("_")
		
		rnum2 = int(tmp2[1])+1 # PROFASI COUNTS RESIDUE NUMBERS FROM 0 !!!
		rname2 = tmp2[2]
		aname2 = tmp2[3].strip("_")
		
		assert aname1 in atomtypes
		assert aname2 in atomtypes
		
		#if potential=="FMULTIGAUSS":
			
		atoms = atomsFromPDB(pdb_files[0])
		dist = distances_list[0][x]
		
		if len(pdb_files)==2:
			atoms2 = atomsFromPDB(pdb_files[1])
			dist2 = distances_list[1][x]
		
		
		
		a1 = None
		a2 = None
		for a in atoms.keys():
			#print atoms[a]
			
			if atoms[a]["rnum"]==rnum1 and atoms[a]["rname"]==rname1 and atoms[a]["aname"]==aname1:
				assert a1==None
				a1 = a
				#print "  found a1",i[0],"in",pdb_files[0]
				nFound += 0.5
			if atoms[a]["rnum"]==rnum2 and atoms[a]["rname"]==rname2 and atoms[a]["aname"]==aname2:
				assert a2==None
				a2 = a
				#print "  found a2",i[1],"in",pdb_files[0]
				nFound += 0.5
		
		
		
		if a1 != None and a2 != None:
			
			aad = atomDist(a1,a2,atoms)
			if len(pdb_files)==1:
				energy = getEnergy(t, parms, aad)
				print pdb_files[0],x,t,aad,energy
				#print "!!!!!!!!!"
				totalEnergy1.append(energy)
			diff = abs(aad-dist)
			if len(pdb_files)==1: 
				if diff >= tolerance: countDiff += 1
				avgAbsDiff.append(diff)
		else:
			if any([xxx!=None for xxx in [a1,a2]]): print a1,a2
		
		if len(pdb_files)==2:
		
			a1 = None
			a2 = None
			for a in atoms2.keys():
				#print atoms[a]
			
				if atoms2[a]["rnum"]==rnum1 and atoms2[a]["rname"]==rname1 and atoms2[a]["aname"]==aname1:
					assert a1==None
					a1 = a
					#print "    found a1",i[0],"in",pdb_files[1]
					nFound2 += 0.5
				if atoms2[a]["rnum"]==rnum2 and atoms2[a]["rname"]==rname2 and atoms2[a]["aname"]==aname2:
					assert a2==None
					a2 = a
					#print "    found a2",i[1],"in",pdb_files[1]
					nFound2 += 0.5
					
			
			
			if a1 != None and a2 != None:
				aad2 = atomDist(a1,a2,atoms2)
				
				energy1 = getEnergy(t, parms,aad)
				energy2 = getEnergy(t, parms,aad2)
				if energy1 > -1.0e-3: energy1 = 0.0
				if energy2 > -1.0e-3: energy2 = 0.0
				
				
				print x,t,aad, aad2,energy1,energy2
				totalEnergy1.append(energy1)
				totalEnergy2.append(energy2)
				diff2 = abs(aad2-dist2)
				
				
				if not(diff < tolerance or diff2 < tolerance): 
					if diff >= tolerance: 
						avgAbsDiff.append(min(diff,diff2))
						countDiff += 1
					
					print diff,diff2
					
	
	print "Comparing %s with XML file:"%(pdb_files[0])
	print "  --> found %d of %i contacts in %s"%(nFound,len(xml_lines),outfilename)
	
	if len(pdb_files)==2:
		print "Comparing %s with XML file:"%(pdb_files[1])
		print "  --> found %d of %i contacts in %s"%(nFound2,len(xml_lines),outfilename)
	
	if len(avgAbsDiff)==0:
		avgAbsDiff = 0
	else:
		avgAbsDiff = sum(avgAbsDiff)/float(len(avgAbsDiff))
	print "  average difference of native distance (%i): %f Angstrom"%(countDiff, avgAbsDiff )
	
	print "done"
	
	
	print "Total energy = ", sum(totalEnergy1), sum(totalEnergy2)
	
	

def getEnergy(potential, parameters, x=None):
	#print parameters
	
	if potential=="LJ":
		assert len(parameters) == 2
		minimum = parameters[1]
		depth = parameters[0]
		if x==None: x = minimum
		
		ratio = minimum / x
		return depth * ( math.pow(ratio, 12.0) - 2 * math.pow(ratio, 6.0 ) ) 
	
	elif potential=="GAUSS":
		assert len(parameters) == 5
		minimum = parameters[0]
		radius = parameters[1]
		steepness = parameters[2]
		width = parameters[3]
		depth = parameters[4]
		if x==None: x = minimum
		
		Rij = steepness * math.pow(radius/x, 12.0)
		
		if width != 0.0:
			Gij = - math.exp( - math.pow(x-minimum, 2.0) / (2 * width * width) )
		else:
			Gij = 0.0
		return  Rij + depth * Gij + Rij * Gij
	
	elif potential=="DUALGAUSS":
		assert len(parameters)==8
		minimum = parameters[0]
		minimum2 = parameters[1]
		radius = parameters[2]
		steepness = parameters[3]
		width = parameters[4]
		depth = parameters[5]
		width2 = parameters[6]
		depth2 = parameters[7]
		
		#print depth,depth2
		
		
		
		Rij = steepness * math.pow(radius/x, 12.0)
		
		if width != 0.0:
			Gij = - math.exp( - math.pow(x-minimum, 2.0) / (2 * width * width) )
			
		else:
			Gij = 0.0
			
		
		if width2 != 0.0:
			Gij2 = - math.exp( - math.pow(x-minimum2, 2.0) / (2 * width2 * width2) )
			
		else:
			Gij2 = 0.0
			
		
		return Rij  +  depth * Gij + Rij * Gij  +  depth2 * Gij2 + Rij * Gij2
	elif potential=="FDUALGAUSS":
		assert len(parameters)==7
		minimum = parameters[0]
		minimum2 = parameters[1]
		radius = parameters[2]
		steepness = parameters[3]
		width = parameters[4]
		depth = parameters[5]
		width2 = parameters[6]
		
		
		#print depth,depth2
		
		
		
		Rij = steepness * math.pow(radius/x, 12.0)
		
		if width != 0.0:
			Gij = - math.exp( - math.pow(x-minimum, 2.0) / (2 * width * width) )
			
		else:
			Gij = 0.0
			
		
		if width2 != 0.0:
			Gij2 = - math.exp( - math.pow(x-minimum2, 2.0) / (2 * width2 * width2) )
			
		else:
			Gij2 = 0.0
			
		
		return depth * ( (1+(1.0/depth)*Rij) * (1+Gij) * (1+Gij2) -1)
	elif potential=="FMULTIGAUSS":
		assert len(parameters)==7
		minima = [float(p) for p in parameters[0].split(",")]
		
		radius = parameters[1]
		steepness = parameters[2]
		width = parameters[3]
		depth = parameters[4]
		
		
		
		#print depth,depth2
		wells = 1
		for m in minima:
			if width != 0.0:
				Gij = - math.exp( - math.pow(x-m, 2.0) / (2 * width * width) )
			
			else:
				Gij = 0.0
			
			wells *= Gij
		
		Rij = steepness * math.pow(radius/x, 12.0)
		
		
			
		
		
			
		
		return depth * ( (1+(1.0/depth)*Rij) * wells -1)


def plotContactMap(n, contacts, inPDBFilename, cutoff, chainSep, atoms, contactType):
	
	import matplotlib.pyplot as p
	import matplotlib.cm as cm
	
	cdict = {}
	for c in contacts:
		a1 = c[0]
		a2 = c[1]
		d = c[2]
		r1 = atoms[a1]["rnum"]
		r2 = atoms[a2]["rnum"]
		cdict[(r1,r2)] = d
	
	matrix = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	
	for i in xrange(n):
		
		for j in xrange(n):
			
			
			
			if (i,j) in cdict:
				matrix[i][j] = cdict[(i,j)]
				
		
	
	title = "%s; cutoff=%s; csep=%i; n=%i\n%s"%(inPDBFilename,str(round(cutoff,2)),chainSep,len(contacts),contactType)

	imax = p.matshow(matrix,cmap=cm.binary,aspect='equal')
	p.title(title)
	cbar = p.colorbar()
	p.xlabel("residue position")
	p.ylabel("residue position")
	cbar.set_label("CA dist")
	#imax.get_axes().set_axisbelow(True)
	p.show()






def getLength(atoms):
	maxrnum = 0
	for a in atoms:
		if atoms[a]["rnum"]>maxrnum:
			maxrnum = atoms[a]["rnum"]
	return maxrnum










if __name__=="__main__":
	##################
	# HANDLE OPTIONS (-h for help)
	description="%prog prepares native contact files in XML format for PROFASI as modified by TS.\n Default values for options are given in parentheses."
	version="2014-09-17; written by Tobias Sikosek"
	usage="\n\tpython %prog [OPTIONS] -s filename.pdb"
	
	op=optparse.OptionParser(description=description, usage=usage, version=version)
	
	op.add_option('-s','--inPDBFilename', action="store", type="string", dest="inPDBFilename",help="structure file in PDB format. Comma-separate multiple files in order to use potential type FMULTIGAUSS and omit -S option.")
	op.add_option('-S','--inPDBFilename2', action="store", type="string", dest="inPDBFilename2",help="second structure file in PDB format")
	op.add_option('-t','--potentialType', action="store", type="string", dest="potentialType",help="type of interaction potential in PROFASI: LJ, GAUSS, DUALGAUSS, FDUALGAUSS, FMULTIGAUSS  (%default)")
	op.add_option('-c','--contactType', action="store", type="string", dest="contactType",help="type of contact definition: CUTOFF, SHADOW, HYBRID, RESIDUE  (%default)")
	op.add_option('-a','--contactAtoms', action="store", type="string", dest="contactAtoms",help="list of atom types for native contacts, e.g. CA,CB,C,N,O  (%default)")
	op.add_option('-n','--chainSep', action="store", type="int", dest="chainSep",help="min number of residues between contacting atoms  (%default)")
	op.add_option('-d','--distCutoff', action="store", type="float", dest="distCutoff",help="distance cutoff for defining contacts (%default Ang)")
	op.add_option('-m','--wellDepth', action="store", type="float", dest="wellDepth",help="well depth at native distance (%default)")
	op.add_option('-w','--wellWidth', action="store", type="float", dest="wellWidth",help="GAUSS/DUALGAUSS: width of well (%default)")
	op.add_option('-M','--wellDepth2', action="store", type="float", dest="wellDepth2",help="well depth at native distance 2 (%default)")
	op.add_option('-W','--wellWidth2', action="store", type="float", dest="wellWidth2",help="GAUSS/DUALGAUSS: width of well 2 (%default)")
	op.add_option('-r','--repulDist', action="store", type="float", dest="repulDist",help="GAUSS/DUALGAUSS: excluded volume radius.  0=no repulsion; -x=radius is min distance per contact - x (FMULTIGAUSS)  (%default)")
	op.add_option('-p','--repulSteep', action="store", type="float", dest="repulSteep",help="GAUSS/DUALGAUSS: steepness of repulsion (%default)")
	op.add_option('-e','--ensembleCutoff', action="store", type="float", dest="ensembleCutoff",help="when PDB file contains multiple models, retain contacts present in this fraction of models (%default)")
	op.add_option('-N','--normEnergy', action="store",type="float", dest="normEnergy",help="normalize contact energies to given amount per structure. Overrides options m and M. (%default=don't normalize)")
	op.add_option('--gmxpath', action="store", type="string", dest="gmxpath",help="path to Gromacs executables (SHADOW)")
	op.add_option('--scmpath', action="store", type="string", dest="scmpath",help="path to SCM.jar  (SHADOW)")
	op.add_option('-b', action='store_true', dest='useBestModel', help='use the model closest to the average contact distances from the ensemble')
	op.add_option('-L','--label',action='store', type="string",dest="label",help="string that will be appended to the keyword 'DistanceRestraint' in the Profasi energy output")
	
	op.set_defaults(inPDBFilename="",
					inPDBFilename2="",
					potentialType="LJ",
					contactType="CUTOFF",
					contactAtoms="CA",
					chainSep=3,
					distCutoff=6.0,
					wellDepth=1.0,
					wellWidth=0.5,
					wellDepth2=1.0,
					wellWidth2=0.5,
					repulDist=0.0,
					repulSteep=1.0,
					ensembleCutoff=0.8,
					normEnergy=0,
					gmxpath="",
					scmpath="",
					useBestModel=False,
					label="")
	
	opt, args = op.parse_args()
								
	
	inPDBFilename = opt.inPDBFilename
	inPDBFilename2 = opt.inPDBFilename2
	potentialType = opt.potentialType
	contactType = opt.contactType
	contactAtoms = opt.contactAtoms
	chainSep = opt.chainSep
	distCutoff = opt.distCutoff
	wellDepth = opt.wellDepth
	wellWidth = opt.wellWidth
	wellDepth2 = opt.wellDepth2
	wellWidth2 = opt.wellWidth2
	repulDist = opt.repulDist
	repulSteep = opt.repulSteep
	ensembleCutoff = opt.ensembleCutoff
	normEnergy = opt.normEnergy
	gmxpath = opt.gmxpath
	scmpath = opt.scmpath
	useBestModel = opt.useBestModel
	label = opt.label
	
	assert potentialType in ["LJ","GAUSS","DUALGAUSS","FDUALGAUSS","FMULTIGAUSS"]
	assert contactType in ["CUTOFF","SHADOW","HYBRID","RESIDUE"]
	assert distCutoff > 1.0, "cutoff in Angstrom!"
	
	
	#if "," in inPDBFilename:
	#	assert inPDBFilename2==""
	#	assert potentialType=="FMULTIGAUSS"
	
	#if potentialType=="FMULTIGAUSS":
	#	assert inPDBFilename2==""
	#	assert "," in inPDBFilename
		
	
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
	
	pdb_files = [inPDBFilename]
	if inPDBFilename2!="": pdb_files.append(inPDBFilename2)
	
	
	
	pdb2contacts = {}
	pdb2atoms = {}
	pdb2models = {}
	
	bestModels = []
	
	atoms_models = []
	
	for p in pdb_files:
		#print pdb_files
		models = splitModels(p)
		print "\n\n%i model(s) found in %s."%(len(models),p)
		models_contact_lists = []
		for m in models:
			contacts = []
			con_c=[]
			con_s=[]
			
			atoms = atomsFromPDB(m)
			
			if contactType=="CUTOFF":
				contacts = getCutoffContacts(pdbfilename=m, atomtypes=contactAtoms, cutoff=distCutoff, chainDist=chainSep, verbose=False, res2SS = None)
			elif contactType=="SHADOW":
				contacts = getShadowContacts(pdbfilename=m, outfilename=p[:-4]+".shadow", cutoff=6.0, chainDist=3, gmxpath=gmxpath, scmpath=scmpath, res2SS = None)
			elif contactType=="HYBRID":
				con_c = getCutoffContacts(pdbfilename=m, atomtypes=contactAtoms, cutoff=distCutoff, chainDist=chainSep, verbose=False, res2SS = None)
				con_s = getShadowContacts(pdbfilename=m, outfilename=p[:-4]+".shadow", cutoff=6.0, chainDist=3, gmxpath=gmxpath, scmpath=scmpath, res2SS = None)
				contacts = mixCutoffShadow(con_c, con_s)
			elif contactType=="RESIDUE":
				contacts = getResidueContacts(pdbfilename=m, cutoff=distCutoff, chainDist=chainSep, atoms=atoms,verbose=False, res2SS = None)
			else:
				print "unknown contact type: ",contactType
				sys.exit(1)
			
			#print "model:",m,"contacts:",len(contacts),"type:",contactType
			models_contact_lists.append(contacts)
			#print len(models_contact_lists[-1])
			
			atoms_models.append(atoms)
		
		pdb2models[p] = models
		
		if potentialType=="FMULTIGAUSS":
			print "FMULTIGAUSS ensemble contacts"
			finalcontacts, ref_atoms, dists, means, stddevs = filterSuperEnsembleContacts(models_contact_lists, atoms_models, contactAtoms, ensembleCutoff, contactType)
			#atoms = atomsFromPDB(models[0])
			#bestModels.append(models[0])
		else:
			finalcontacts, bestModelIndex = filterEnsembleContacts(models_contact_lists, atoms, contactAtoms, ensembleCutoff, contactType)
		
		if useBestModel:
			atoms = atomsFromPDB(models[bestModelIndex])
			bestModels.append(models[bestModelIndex])
			consStatement += "_m%i"%(bestModelIndex+1)
		else:
			atoms = atomsFromPDB(models[0])
			bestModels.append(models[0])
		
		if potentialType=="FMULTIGAUSS":
			pdb2atoms[p] = atoms_models
			#pdb2contacts[p] = models_contact_lists
		else:
			pdb2atoms[p] = atoms
		
		pdb2contacts[p] = contacts
		
		
	
	pbond_legend = None
	
	if len(pdb_files)==2:
		if contactType=="LJ": removeShared = True
		else: removeShared = False
		pbond_legend = []
		contacts = mixContactsBistable(pdb2contacts[pdb_files[0]], pdb2contacts[pdb_files[1]], removeShared=removeShared, legend=pbond_legend)
	
	outatoms = "".join([i for i in contactAtoms])
	
	if normEnergy==0:
		if len(pdb_files) == 1:
			totalEnergy = wellDepth*len(pdb2contacts[inPDBFilename])
		elif len(pdb_files) == 2:
			totalEnergy = wellDepth*len(pdb2contacts[inPDBFilename]) + wellDepth2*len(pdb2contacts[inPDBFilename2])
		print "Total energy (approx.):",totalEnergy
	else:
		totalEnergy = normEnergy * len(pdb_files)
		print "Total energy:",totalEnergy
	
	if potentialType in ["GAUSS","FMULTIGAUSS"]:
		width_label = "_w"+str(round(wellWidth,1))
		repul_label = "_R"+str(round(repulDist,1))
	else:
		width_label = ""
		repul_label = ""
	
	
	
	outfilename = "%s_%s_%s_%s_e%s_c%s_n%s%s%s_%s%s_nativeContacts.xml"%(inPDBFilename.split("/")[-1],inPDBFilename2.split("/")[-1],potentialType,contactType,totalEnergy,str(distCutoff),str(chainSep),width_label,repul_label,outatoms,consStatement)
	
	if contactType=="HYBRID" and not 'CA' in contactAtoms:
		contactAtoms.append('CA')
	
	
	
	if potentialType=="LJ": 		 writeProfasiContactsLJ(finalcontacts, pdb2atoms[inPDBFilename], outfilename, wellDepth, norm=normEnergy,label=label)
	elif potentialType=="GAUSS": 	 writeProfasiContactsGAUSS(finalcontacts, pdb2atoms[inPDBFilename], outfilename, repulDist, repulSteep, wellWidth, wellDepth, norm=normEnergy,label=label)
	elif potentialType=="DUALGAUSS": writeProfasiContactsDUALGAUSS(pdb2contacts[pdb_files[0]], pdb2contacts[pdb_files[1]], pdb2atoms[inPDBFilename],pdb2atoms[inPDBFilename2], contactAtoms, outfilename, repulDist, repulSteep, wellWidth, wellDepth, wellWidth2, wellDepth2, norm=normEnergy,label=label)
	elif potentialType=="FDUALGAUSS": writeProfasiContactsFDUALGAUSS(pdb2contacts[pdb_files[0]], pdb2contacts[pdb_files[1]], pdb2atoms[inPDBFilename],pdb2atoms[inPDBFilename2], contactAtoms, outfilename, repulDist, repulSteep, wellWidth, wellDepth, wellWidth2,label=label)
	elif potentialType=="FMULTIGAUSS": writeProfasiContactsFMULTIGAUSS(finalcontacts, dists, ref_atoms, outfilename, repulDist, repulSteep, wellWidth, wellDepth, norm=normEnergy,label=label)
	
	if potentialType=="FMULTIGAUSS":
	#	for m in xrange(len(models)):
		writeChimeraPseudobonds(finalcontacts, ref_atoms, outfilename[:-4]+".pseudobonds",legend=pbond_legend)
	#else:
	writeChimeraPseudobonds(finalcontacts, atoms, outfilename[:-4]+".pseudobonds",legend=pbond_legend)
	
	print
	if not potentialType=="FMULTIGAUSS":
		verifyProfasiXML(outfilename, bestModels, contactAtoms)
	
	if False:
		nres = getLength(atoms)
		plotContactMap(nres, finalcontacts, inPDBFilename, distCutoff, chainSep, atoms, contactType)
	
	
	
	