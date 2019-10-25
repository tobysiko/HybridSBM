#!/usr/bin/python
import sys
import numpy
import optparse
import math
import subprocess
import os, os.path
import copy
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import operator

if True: #some lists
	BB_atoms = ["C","CA","O","N","OXT","CB"]
	AA_3_to_1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
	AA_1_to_3 = {'G':'GLY','A':'ALA','V':'VAL','L':'LEU','I':'ILE','C':'CYS','M':'MET','F':'PHE','Y':'TYR','W':'TRP','P':'PRO','S':'SER','T':'THR','N':'ASN','Q':'GLN','D':'ASP','E':'GLU','H':'HIS','K':'LYS','R':'ARG'}
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

def splitModels(pdbfilename, labels=[]):
	models = []
	model_nr = 0
	atomlines = []
	
	for line in open(pdbfilename).readlines():
		linetype = line[:6]
		
		if linetype == "MODEL ":
			model_nr = int(line.replace("MODEL",""))
			#print "model nr ",model_nr
			atomlines = []
		elif "REMARK original_filename" in line:
			labels.append(os.path.splitext(line.strip().split()[2])[0])
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
	
	

def filterSuperEnsembleContacts(models_contact_lists, models_atoms, model_chains, atypes, ensembleCutoff, contactType, verbose=False, map2msa=None, modelSeqMsaIDs=[], isDimer=False):
	if verbosity >= 1: print "super ensemble contacts"
	contactsDict = {}
	
	atomtypes = copy.copy(atypes)
	
	if contactType == "HYBRID" and not "CA" in atomtypes:
		atomtypes  += ['CA']
	
	nModels = len(models_contact_lists)
	contactsPerModel = [[] for i in xrange(nModels)]
	
	for m in xrange(nModels):
		model_contacts = models_contact_lists[m]
		if verbosity >= 2: print "%i contacts in model %s"%(len(model_contacts),m)
		
		modelatoms = models_atoms[m]
		if verbosity >= 2: print len(modelatoms)," atoms in modelatoms"
		
		if map2msa != None and verbosity >= 2: print "NOTE: contact positions refer to MSA!!!"
		
		for c in model_contacts:
			a1 = c[0]
			a2 = c[1]
			d  = c[2]
			if len(c) == 5:
				ch1 = c[3]
				ch2 = c[4]
			else:
				ch1 = None
				ch2 = None
			
			#print c
			assert modelatoms[a1]["aname"] in atypes, modelatoms[a1]
			assert modelatoms[a2]["aname"] in atypes, modelatoms[a2]
			
			atype1 = modelatoms[a1]["aname"]
			atype2 = modelatoms[a2]["aname"]
			
			if map2msa != None: # map positions to MSA
				assert modelSeqMsaIDs[m] < len(map2msa)
				msalength = max( map2msa[ modelSeqMsaIDs[m] ])
				
				if modelatoms[a1]["rnum"]>msalength or modelatoms[a2]["rnum"]>msalength:
					print "WARNING: contact position exceeds MSA length! Skip.",(modelatoms[a1]["rnum"],modelatoms[a2]["rnum"])
					continue
				else:
					assert modelatoms[a1]["rnum"] in map2msa[ modelSeqMsaIDs[m] ], (modelatoms[a1]["rnum"], map2msa[ modelSeqMsaIDs[m] ])
					assert modelatoms[a2]["rnum"] in map2msa[ modelSeqMsaIDs[m] ], (modelatoms[a2]["rnum"], map2msa[ modelSeqMsaIDs[m] ])
				rnum1 = map2msa[ modelSeqMsaIDs[m] ][ modelatoms[a1]["rnum"] ]
				rnum2 = map2msa[ modelSeqMsaIDs[m] ][ modelatoms[a2]["rnum"] ]
				if verbosity >=2: print "Mapped contact (%i,%i) in sequence %s to (%i,%i) in alignment."%(modelatoms[a1]["rnum"], modelatoms[a2]["rnum"], modelSeqMsaIDs[m], rnum1, rnum2)
				if not isDimer:
					mprcd = minPairwiseResidueCADist(modelatoms, modelatoms[a1]["rnum"], modelatoms[a2]["rnum"])
					
				else:
					mprcd = atomDist(a1,a2,modelatoms)
				assert d == mprcd, str(d)+" != "+str(mprcd)
			else:
				rnum1 = modelatoms[a1]["rnum"]
				rnum2 = modelatoms[a2]["rnum"]
			
			if ch1 != None and ch2 != None:
				# write contacts to contact dictionary
				if (rnum1,rnum2,atype1,atype2,ch1,ch2) in contactsDict:
				
					if len(contactsDict[(rnum1,rnum2,atype1,atype2,ch1,ch2)]) < nModels:
						contactsDict[(rnum1,rnum2,atype1,atype2,ch1,ch2)].append(d)
					#print (rnum1,rnum2,atype1,atype2)
				
				else:
					contactsDict[(rnum1,rnum2,atype1,atype2,ch1,ch2)] = [d]
					#print "new:", (rnum1,rnum2,atype1,atype2)
			else:	
				# write contacts to contact dictionary
				if (rnum1,rnum2,atype1,atype2) in contactsDict:
				
					if len(contactsDict[(rnum1,rnum2,atype1,atype2)]) < nModels:
						contactsDict[(rnum1,rnum2,atype1,atype2)].append(d)
					#print (rnum1,rnum2,atype1,atype2)
				
				else:
					contactsDict[(rnum1,rnum2,atype1,atype2)] = [d]
					#print "new:", (rnum1,rnum2,atype1,atype2)
	
	# which contacts are found in all models? what are the means and standard deviations of distances?
	
	#print len(models_contact_lists),len(model_contacts),len(contactsDict)
	
	counter = 0
	finalContacts = []
	dists = []
	global_min_dist = 10000
	global_max_dist = -1
	min_range = 100000
	max_range = -1
	
	if verbosity >= 1: print "Collected %i contacts from all models."%len(contactsDict.keys())
	
	fracHistoSize = 100
	fracHisto = [0 for h in xrange(fracHistoSize)]
	
	for i in sorted(contactsDict.keys()):
		#print i,contactsDict[i]
		contactCount = len(contactsDict[i]) 
		
		assert contactCount <= nModels, str([contactCount, nModels])
		
		
		
		frac = contactCount / float(nModels)
		
		for h in xrange(int(frac*100)):
			fracHisto[h] +=1
		
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
			
			if verbosity >=2: print "(%s) accepted atom pair:"%str(round(frac,1)), i, " - ", contactCount, "models with this contact.", "Min,Max:",min_d,max_d," Range:",range_d
			
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
				if verbosity >= 1: print "atom is not allowed for contact!",ref_atoms[i[0]]["aname"],ref_atoms[i[1]]["aname"]
		else:
			if verbosity >=2: print "(%s) REJECTED atom pair:"%str(round(frac,1)), i, " - ", contactCount, "models with this contact."
	print
	print "Consensus fraction contact histogram:"
	print "fraction of models","\t","n contacts"
	for h in xrange(fracHistoSize):
		print (h+1.0)/fracHistoSize,"\t", fracHisto[h]
		
	print
	
	if verbosity >= 1: 
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

def getSequence(atoms):
	rd={}
	seq=""
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in rd:
			rd[rnum].append(atoms[a])
		else:
			rd[rnum]=[atoms[a]]
			seq += AA_3_to_1[atoms[a]["rname"].upper()]
	
	return seq
def getChainSequences(atoms):
	rd={}
	seqs=[] # list of sequences, one per chain
	seq=""
	currchain = ""
	for a in atoms:
		rnum = atoms[a]["rnum"]
		ch = atoms[a]["chain"]
		if ch != currchain:
			rd={}
			if seq != "": 
				seqs.append(seq)
				seq = ""
				currchain = ch
			else:
				currchain = ch
		if rnum in rd:
			rd[rnum].append(atoms[a])
		else:
			rd[rnum]=[atoms[a]]
			seq += AA_3_to_1[atoms[a]["rname"].upper()]
	if seq != "":
		seqs.append(seq)
	return seqs

def getResidueDict(atoms):
	rd = {}
	
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in rd:
			rd[rnum].append(atoms[a])
		else:
			rd[rnum]=[atoms[a]]
	
	
	return rd
def getChainResidueDict(atoms):
	chains ={}
	
	for a in atoms:
		ch = atoms[a]["chain"]
		rnum = atoms[a]["rnum"]
		if ch in chains:
			if rnum in chains[ch]:
				chains[ch][rnum].append(atoms[a])
			else:
				chains[ch][rnum]=[atoms[a]]
		else:
			chains[ch] = {}
			chains[ch][rnum]=[atoms[a]]
	
	return chains
def minPairwiseResidueDist(atoms, res1,res2):
	residues = getResidueDict(atoms)
	
	mindist = 10e10
	
	for ai in residues[res1]:
		if 'H' in ai["aname"]: continue
		for aj in residues[res2]:
			if 'H' in aj["aname"]: continue
			d = atomDist(ai["anum"],aj["anum"],atoms)
			if d<mindist:
				mindist=d
	return mindist
def minPairwiseResidueCADist(atoms, res1, res2):
	residues = getResidueDict(atoms)
	for ai in residues[res1]:
		for aj in residues[res2]:
			if ai["aname"]=="CA" and aj["aname"]=="CA":
				return atomDist(ai["anum"],aj["anum"],atoms)
			
	return None
def getInterchainResidueContacts(atoms, chainsep, cutoff, verbose=False, selection=[], selectExclusive=False, interChainContacts=True, intraChainContacts=True):

	contacts = []
	chainresidues = getChainResidueDict(atoms)
	nchains = len(chainresidues)
	if verbosity >= 2: 
		print nchains,"chains (getResidueContacts)"
		for i in xrange(nchains):
			print "%s: %i residues"%(chainresidues.keys()[i], len(chainresidues[chainresidues.keys()[i]]))
	
	
	for ch1 in xrange(nchains):
		ch1key = chainresidues.keys()[ch1]
		chlen1 = len(chainresidues[ch1key])
		
		for ch2 in xrange(nchains):
			ch2key = chainresidues.keys()[ch2]
			chlen2 = len(chainresidues[ch2key])
			if ch1<=ch2:
				for i in xrange(chlen1):
					if selectExclusive and not i+1 in selection: continue
					if not i+1 in chainresidues[ch1key]: 
						print i+1
						continue
					ri = chainresidues[ch1key][i+1]
					for j in xrange(chlen2):
						
						if ch1!=ch2 or (ch1==ch2 and i+chainsep<j):
							
							if selection!=[] and selectExclusive and not j+1 in selection: continue
		
							if selection!=[] and not selectExclusive and not i+1 in selection and not j+1 in selection: continue
		
							if not j+1 in chainresidues[ch2key]: 
								print j+1
								continue
							rj = chainresidues[ch2key][j+1]
	
					
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
								a1 = getAtomIDforChainResidueContact(i+1,'CA',atoms,ch1key)
								a2 = getAtomIDforChainResidueContact(j+1,'CA',atoms,ch2key)
								
								if (interChainContacts and ch1key != ch2key) or (intraChainContacts and ch1key == ch2key):
									contacts.append( [a1,a2,cadist,ch1key,ch2key] )
								#print i+1,j+1,[a1,a2,cadist,ch1key,ch2key],chlen1,chlen2
								#print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	
	
	return contacts
def getChainIDs(atoms):
	ch={}
	for a in atoms:
		if not atoms[a]["chain"] in ch:
			ch[atoms[a]["chain"]] = True
	return ch.keys()
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
	if verbosity >= 2: print n,"residues (getResidueContacts)"
	
	for i in xrange(n):
		if selectExclusive and not i+1 in selection: continue
		if not i+1 in residues: continue
		ri = residues[i+1]
		for j in xrange(n):
			if selection!=[] and selectExclusive and not j+1 in selection: continue
			
			if selection!=[] and not selectExclusive and not i+1 in selection and not j+1 in selection: continue
			
			if not j+1 in residues: continue
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
					
					contacts.append( [a1,a2,cadist,ch1key,ch2key] )
					#print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	
	
	return contacts
def getAtomIDforResidueContact(resid,atype,atoms):
	
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype:
			
			return a
	
	
	return None

def getAtomIDforChainResidueContact(resid,atype,atoms,chkey):
	
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype and atoms[a]["chain"] == chkey:
			
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
	if verbosity >= 2: print "Running command for SHADOW CONTACTS: "+scm
	sub = subprocess.Popen(scm,stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	contacts = []
	
	if verbosity != 2: 
		print out
		print err
	
	gro2pdb = gro2pdbAtomMap(grofilename, pdbfilename, atomtypes=['CA'])
	
	if not os.path.exists(outfilename):
		print "SHADOW contacts were not written! Have you checked that 'SCM.jar' is available? Use '--scmpath' option."
		sys.exit(1)
	
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

def getResNameFromNumber(atoms, rnum, aname, chain=None):
	for a in atoms:
		
		#print atoms[a]["rnum"], rnum, atoms[a]["aname"] , aname, atoms[a]["chain"], chain
		if (chain!=None and atoms[a]["rnum"]==rnum and atoms[a]["aname"] ==aname and atoms[a]["chain"]==chain) or \
		 (chain==None and atoms[a]["rnum"]==rnum and atoms[a]["aname"] ==aname):
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

def writeProfasiContactsFMULTIGAUSS(contacts, atoms, atomtypes, outfilename, radius, steepness, width, depth, norm=0,label="",includeNonNative=False, chainSep=3, smoothGaussian=False, restraintType="distance_restraints", mapMSA2Seq=None, nchains=1, chainmap={"A":0, "B":1, "C":2, "D":3, "E":4, "F":5}, generateChainInfo = False, aligned_sequence=None,genericSequenceOfLength=None):
	
	if aligned_sequence != None:
		sequence = degap(aligned_sequence)
	elif genericSequenceOfLength != None:
		sequence = "G"*genericSequenceOfLength
	else:
		sequence = getSequence(atoms)
	#print mapMSA2Seq
	
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
 <data>\n"""%(restraintType,label))
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
 <data>\n"""%(restraintType,label))
	
	seqdists = 0
	#print radius, steepness, width, depth
	
	nskipped = 0
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		
		msa_rnum1 = c[0]
		msa_rnum2 = c[1]
		
		if mapMSA2Seq != None:
			
			
			if msa_rnum1 in mapMSA2Seq:
				rnum1 = mapMSA2Seq[ msa_rnum1 ]
			else:
				if verbosity >= 2: print "Could not map contact %s to this sequence, because aligned position %i has no match in this sequence."%(str(c),c[0])
				nskipped += 1
				continue
			
			if msa_rnum2 in mapMSA2Seq:
				rnum2 = mapMSA2Seq[ msa_rnum2 ]
			else:
				if verbosity >= 2: print "Could not map contact %s to this sequence, because aligned position %i has no match in this sequence."%(str(c),c[1])
				nskipped += 1
				continue
			if verbosity >= 2: print "Mapped aligned contact positions (%i,%i) to sequence positions (%i,%i)."%(c[0],c[1],rnum1,rnum2)
		else:
			
			rnum1 = c[0]
			rnum2 = c[1]
		
		
		
		#print c,rnum1,rnum2
		atype1 = c[2]
		atype2 = c[3]
		
		
		if len(c) == 7:
			
			ch1 = c[4]
			ch2 = c[5]
			dists = c[6]
		else:
			dists = c[4]
			ch1 = None
			ch2 = None
		
		energies = [getEnergy(i, dists, radius, steepness, width, depth) for i in dists]
		assert all([i==-depth for i in energies])
		
		if False: # FOR DEBUGGING
			residues=getResidueDict(atoms)
			if not rnum1 in residues:
				print "Could not find",rnum1
				continue
			if not rnum2 in residues:
				print "Could not find",rnum2
				continue
			a1 = getAtomIDforResidueContact(rnum1, atype1, atoms)
			a2 = getAtomIDforResidueContact(rnum2, atype2, atoms)
			d1 = round(minPairwiseResidueDist(atoms, rnum1, rnum2),2)
			#d2 = round(atomDist(a1,a2,verifyPDBatoms),2)
			d3 = round(minPairwiseResidueCADist(atoms, rnum1, rnum2),2)
			#print verifyPDBatoms[a1]
			#print verifyPDBatoms[a2]
			#print dists,d3,"---",d1,distCutoff
			#print d1,d2,d3
			#assert d3 >= min(dists) and d3 <= min(dists)
			#if d1 > distCutoff: 
			#	print "Residue distance above cutoff",d1,distCutoff
			#	continue
			if d3 >= min(dists) and d3 <= min(dists):
				print "CA distance outside range",d3,dists
				print
		
		mindist = min(dists)
		if radius < 0:
			final_radius = mindist + radius
		else:
			final_radius = radius
		#print rnum1,rnum2,atype1,atype2
		
		#rname1 = getResNameFromNumber(atoms, msa_rnum1, atype1, chain=ch1)
		#rname2 = getResNameFromNumber(atoms, msa_rnum2, atype2, chain=ch2)
		
		#assert rname1 != None and rname2 != None, (c, rname1,rname2)
		
		#if rname1 == None or rname2 == None:
			# in the MSA there is no amino acid at this position, but only a gap '-'
			#nskipped += 1
			#continue
		#	pass
		#else:
		
		if aligned_sequence != None: assert aligned_sequence[msa_rnum1-1] !='-' and aligned_sequence[msa_rnum2-1] !='-'
		
		#if aligned_sequence[msa_rnum1-1] =='-' or aligned_sequence[msa_rnum2-1] =='-':
		#	nskipped += 1
		#	print "HOW?!"
		#	continue
		##rnum1 = mapMSA2Seq[rnum1]
		#rnum2 = mapMSA2Seq[rnum2]
		rname1 = AA_1_to_3[sequence[rnum1-1]].upper()
		rname2 = AA_1_to_3[sequence[rnum2-1]].upper()
		#print rnum1,mapMSA2Seq[rnum1], rname1
		#print rnum2,mapMSA2Seq[rnum2], rname2
	
		seqdists += abs(rnum1-rnum2)
		
		if (ch1 == None and ch2 == None) or generateChainInfo:
		
			for chain in xrange(nchains):
				atomstring1 = "%i/%s/%s/_%s_"%(chain, str(rnum1-1), rname1, atype1)#profasiAtomString(a1, atoms)
				atomstring2 = "%i/%s/%s/_%s_"%(chain, str(rnum2-1), rname2, atype2)#profasiAtomString(a2, atoms)
		
				if smoothGaussian:
					out.write( "  %s  %s  FMULTIGSMOOTH  %s  %s  %s %s %s %s\n"%(atomstring1, atomstring2, str(min(dists)), str(max(dists)), str(final_radius), str(steepness), str(width), str(depth) ) )
				else:	
					out.write( "  %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n"%(atomstring1, atomstring2, ",".join([str(i) for i in dists]), str(final_radius), str(steepness), str(width), str(depth) ) )
		else:
			chid1 = chainmap[ch1]
			chid2 = chainmap[ch2]
			atomstring1 = "%i/%s/%s/_%s_"%(chid1, str(rnum1-1), rname1, atype1)#profasiAtomString(a1, atoms)
			atomstring2 = "%i/%s/%s/_%s_"%(chid2, str(rnum2-1), rname2, atype2)#profasiAtomString(a2, atoms)
	
			if smoothGaussian:
				out.write( "  %s  %s  FMULTIGSMOOTH  %s  %s  %s %s %s %s\n"%(atomstring1, atomstring2, str(min(dists)), str(max(dists)), str(final_radius), str(steepness), str(width), str(depth) ) )
			else:	
				out.write( "  %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n"%(atomstring1, atomstring2, ",".join([str(i) for i in dists]), str(final_radius), str(steepness), str(width), str(depth) ) )
			
	
	N=getLength(atoms)
	co = seqdists/float(n*N)
	if verbosity >= 2: print "CO,distSum,nContacts,nRes = ", seqdists/float(n*N), seqdists, n, N
	
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
 </%s>"""%(restraintType))
	
	if verbosity >= 2: print "skipped %i out of %i contacts!"%(nskipped, len(contacts))
	
	return outfilename, co


def writeChimeraPseudobonds(contacts, outfilename, legend=None, mapMSA2Seq=None, verifyPDBatoms=None, nchains=1, chainmap={0:"A", 1:"B",2:"C",3:"D",4:"E",5:"F"}, generateChainInfo = False):
	outfile = open(outfilename, 'w')
	
	if legend != None: assert len(contacts) == len(legend),str(len(contacts))+" "+str(len(legend))
	
	maxrange = 0
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		if len(c)==7:
			dists = c[6]
		else:
			dists = c[4]
			
		rn = max(dists) - min(dists)
		if rn > maxrange:
			maxrange = rn
	if verbosity >= 2: print "max.range of CA distances as label in pseudobond file:",maxrange
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		if mapMSA2Seq != None:
			if c[0] in mapMSA2Seq:
				rnum1 = mapMSA2Seq[ c[0] ]
			else:
				#if verbosity >=2: print "Could not map contact %s to this sequence, because aligned position %i has no match in this sequence."%(str(c),c[0])
				continue
			
			if c[1] in mapMSA2Seq:
				rnum2 = mapMSA2Seq[ c[1] ]
			else:
				#if verbosity >=2: print "Could not map contact %s to this sequence, because aligned position %i has no match in this sequence."%(str(c),c[1])
				continue
			#if verbosity >=2: print "Mapped aligned contact positions (%i,%i) to sequence positions (%i,%i)."%(c[0],c[1],rnum1,rnum2)
		else:
			rnum1 = c[0]
			rnum2 = c[1]
		atype1 = c[2]
		atype2 = c[3]
		
		if len(c)==7:
			ch1 = c[4]
			ch2 = c[5]
			dists = c[6]
		else:
			dists = c[4]
			ch1 = None
			ch2 = None
		
		if False:# verifyPDBatoms!=None: # FOR DEBUGGING
			residues=getResidueDict(verifyPDBatoms)
			if not rnum1 in residues:
				print "Could not find",rnum1
				continue
			if not rnum2 in residues:
				print "Could not find",rnum2
				continue
			a1 = getAtomIDforResidueContact(rnum1, atype1, verifyPDBatoms)
			a2 = getAtomIDforResidueContact(rnum2, atype2, verifyPDBatoms)
			d1 = round(minPairwiseResidueDist(verifyPDBatoms, rnum1, rnum2),2)
			#d2 = round(atomDist(a1,a2,verifyPDBatoms),2)
			d3 = round(minPairwiseResidueCADist(verifyPDBatoms, rnum1, rnum2),2)
			#print verifyPDBatoms[a1]
			#print verifyPDBatoms[a2]
			#print dists,d3,"---",d1,distCutoff
			#print d1,d2,d3
			#assert d3 >= min(dists) and d3 <= min(dists)
			if d1 > distCutoff: 
				print "Residue distance above cutoff",d1,distCutoff
				continue
			if d3 >= min(dists) and d3 <= min(dists):
				print "CA distance outside range",d3,dists
				print
		
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
		
		if ch1 == None and ch2 == None:
			pseudobond = "#0:%s@%s #0:%s@%s %s %s"%(str(rnum1), atype1.lower(), str(rnum2), atype2.lower(), color, label)
			outfile.write(pseudobond + "\n")
		else:
			if generateChainInfo:
				for i in xrange(nchains):
					ch1 = chainmap[i]
					ch2 = chainmap[i]
					pseudobond = "#0:%s.%s@%s #0:%s.%s@%s %s %s"%(str(rnum1),ch1, atype1.lower(), str(rnum2),ch2, atype2.lower(), color, label)
					outfile.write(pseudobond + "\n")
			else:
				pseudobond = "#0:%s.%s@%s #0:%s.%s@%s %s %s"%(str(rnum1),ch1, atype1.lower(), str(rnum2),ch2, atype2.lower(), color, label)
				outfile.write(pseudobond + "\n")
				
		

	outfile.close()



def sequencePosMarker(length):
	r = len(str(length))
	
	rows = [["0" for j in xrange(length)] for i in xrange(r)]
	
	for i in range(0, length):
		rev = str(i + 1)[::-1]
		for c in xrange(len(rev)):
			char = rev[c]
			rows[c][i] = char
	
	return "\n".join(["".join(i) for i in rows][::-1])
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


def plotContactMap(contacts, inPDBFilename, cutoff, chainSep, atoms, contactType, fname="contactmap.png",dpi=100, useMSAPos=None):
	if useMSAPos==None:
		n = getLength(atoms)
	else:
		n = useMSAPos
	
	cdict = {}
	for c in contacts:
		if len(c) ==5:
			r1 = c[0]
			r2 = c[1]
			dists = c[4]
		else:
			r1 = c[0]
			r2 = c[1]
			dists = c[6]
		
		#cdict[(r1,r2)] = numpy.mean(dists) # plot MEAN CA distances
		if useMSAPos == None:
			if len(dists)==1:
				cdict[(r1,r2)] = 1.0
			else:
				cdict[(r1,r2)] = max(dists)-min(dists)  # plot RANGE of CA distances
		else:
			cdict[(r1,r2)] = len(dists)
	
	matrix = [[0.0 for j in xrange(n+1)] for i in xrange(n+1) ]
	
	cdict_max = max(cdict.iterkeys(), key=(lambda key:cdict[key]))#max(cdict.iteritems(), key=operator.itemgetter(1))[0]
	#cmap = cm.binary
	for i in xrange(n+1):
		
		for j in xrange(n+1):
			
			if (i,j) in cdict:
				matrix[i][j] = cdict[(i,j)]
				matrix[j][i] = cdict[cdict_max]#cdict[(i,j)]
				#p.scatter([i+1],[j+1],marker="s",s=5,cmap=cmap)
	
	open(fname+".data","w").write("\n".join( [",".join([str(j) for j in i[1:]]) for i in matrix[1:]]))
	
	matrix = [[matrix[i][j] if matrix[i][j] != "NaN" else 0.0 for j in xrange(n+1)] for i in xrange(n+1) ]
	
	title = "%s; cutoff=%s; csep=%i; n=%i\n%s"%(inPDBFilename,str(round(cutoff,2)),chainSep,len(contacts),contactType)

	imax = p.matshow(matrix,cmap=cm.binary,aspect='equal', origin='lower')
	#cmap=cm.afmhot
	
	p.plot([1,n],[1,n],color='k')
	p.title(title)
	cbar = p.colorbar()
	p.xlabel("residue position")
	p.ylabel("residue position")
	p.xlim(0.5,n+0.5)
	p.ylim(0.5,n+0.5)
	p.grid(color="0.5", linestyle=':', linewidth=1.0)
	#cbar.set_label("CA dist")
	if useMSAPos == None:
		cbar.set_label("CA dist. range")
	else:
		cbar.set_label("n models")
	#imax.get_axes().set_axisbelow(True)
	p.savefig(fname,dpi=dpi)
	
	#p.show()
	if verbosity >=2: print "wrote contact map PNG:",fname
	p.close()






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
def plotMultiGaussianContacts(contacts, eps, w, fname, smoothGaussian=False, dpi=100):
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
		
		if len(c) == 7:
			ch1 = c[4]
			ch2 = c[5]
			dists = c[6]
		else:
			dists = c[4]
			ch1 = None
			ch2 = None
		if min(dists) < globmin: globmin = min(dists)
		if max(dists) > globmax: globmax = max(dists)
		
		ax = p.subplot(sq, sq, ci+1)
		#ax.set_title("%i: (%i,%i)"%(ci+1,r1,r2),size=4)
		
		x = numpy.arange(math.floor(globmin)-1.0,math.ceil(globmax)+1.0,0.01)
		ax.set_ylim(-eps-0.1,0.2)
		ax.set_xlim(math.floor(globmin)-1.0,math.ceil(globmax)+1.0)
		
		
		
		if ci in [i for i in xrange(len(contacts)) if i==0 or i%sq==0 ]: 
			ax.set_yticklabels([0.2,0.0,-0.2,-0.4,-0.6,-0.8,-1.0],size=4)
			#ax.set_yticklabels([i.get_text() for i in ax.get_yticklabels()],size=4)
			ax.set_ylabel("E_ij",size=4)
		else: 
			#ax.set_yticks([])
			ax.set_yticklabels([])
		
		if ci in [i for i in range(len(contacts)-sq, len(contacts) ) ]:	
			ax.set_xticklabels(range(int(math.floor(globmin)-1.0),int(math.ceil(globmax)+1.0),1), size=4)
			#ax.set_xticklabels([i.get_text() for i in ax.get_xticklabels()],size=4)
			ax.set_xlabel("r_ij",size=4)
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
		if ch1 == None and ch2 == None:
			p.text(3, 0.3, "%i: (%i,%i)"%(ci+1,r1,r2),size=4)
		else:
			p.text(3, 0.3, "%i: (%i.%s,%i.%s)"%(ci+1,r1,ch1,r2,ch2),size=4)
	
	p.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.4)
	p.savefig(fname,dpi=dpi)
	#p.show()
	if verbosity >=2: print "wrote potentials PNG:",fname
	p.close()

def parseFasta(f):
	ids=[]
	seqs=[]
	l = 0
	need_aln=False
	for seq_record in SeqIO.parse(f, "fasta"):
		ids.append(seq_record.id)
		s=str(seq_record.seq)
		seqs.append(s)
		if l==0: 
			l=len(s)
		else:
			if l != len(s): need_aln = True
	return ids,seqs,need_aln

def mapSeq2MSA(withgaps): # all sequence positions start at 1
	if verbosity >= 2: print "ungapped -> aligned"
	nogaps = withgaps.replace("-","")
	d={}
	assert len(withgaps) >= len(nogaps)
	nogap_counter = 1
	for i in xrange(len(withgaps)):
		if not withgaps[i]=="-":
			d[nogap_counter]=i+1
			if verbosity >= 2: print  nogap_counter,nogaps[nogap_counter-1],i+1,withgaps[i+1-1]
			nogap_counter += 1
	return d
def mapMSA2Seq(withgaps):
	if verbosity >= 2: print "aligned -> ungapped"
	nogaps = withgaps.replace("-","")
	d={}
	assert len(withgaps) >= len(nogaps)
	nogap_counter = 1
	for i in xrange(len(withgaps)):
		if not withgaps[i]=="-":
			d[i+1]=nogap_counter
			if verbosity >= 2: print i+1,withgaps[i+1-1], nogap_counter,nogaps[nogap_counter-1]
			nogap_counter += 1
	return d
	
def mapMSASeq1Seq2(msaseq1,msaseq2):
	if verbosity >= 2: 
		print msaseq1
		print msaseq2
	assert len(msaseq1) == len(msaseq2)
	d={}
	counter1 = 1
	counter2 = 1
	for i in xrange(len(msaseq1)):
		if verbosity >= 2: print i,
		if msaseq1[i] != '-' and msaseq2[i] != '-':
			d[counter1] = counter2
			if verbosity >= 2: print counter1,counter2
			counter1 += 1
			counter2 += 1
		else:
			if verbosity >= 2: print ""
		if msaseq1[i] == '-': counter1 += 1
		if msaseq2[i] == '-': counter2 += 1
	if verbosity >= 2: print "Aligned positions:",len(d)
	return d
def mapMSASeq1ToDegappedSeq2(msaseq1,msaseq2):
	if verbosity >= 2: 
		print msaseq1
		print msaseq2
	
	msa2seq = mapMSA2Seq(msaseq2)
	
	assert len(msaseq1) == len(msaseq2)
	d={}
	counter1 = 1
	counter2 = 1
	for i in xrange(len(msaseq1)):
		if verbosity >= 2: print i,
		if msaseq1[i] != '-' and msaseq2[i] != '-':
			if counter2 in msa2seq:
				d[counter1] = msa2seq[counter2]
			if verbosity >= 2: print counter1,counter2
			counter1 += 1
			counter2 += 1
		else:
			if verbosity >= 2: print ""
		if msaseq1[i] == '-': counter1 += 1
		if msaseq2[i] == '-': counter2 += 1
	if verbosity >= 2: print "Aligned positions:",len(d)
	return d
def reorderFastaByLabel(filename, original_ids, wrong_ids, wrong_ids_seqs):
	newfile = open(filename,"w")
	for o in original_ids:
		index = wrong_ids.index(o)
		newfile.write(">%s\n%s\n"%(wrong_ids[index], wrong_ids_seqs[index]))
	newfile.close()
	
def alignSequencesAndMap(msafilename):
	# if all sequences have same length, no alignment necessary
	o_ids, o_seqs, need_aln = parseFasta(msafilename)
	
	if need_aln:
		alnfile = os.path.splitext(msafilename)[0]+".aln"
		cline = MuscleCommandline("muscle", input=msafilename, out=alnfile)
		stdout, stderr = cline()
		if verbosity >=2: 
			print stdout
			print stderr
		ids,seqs,need_aln = parseFasta(alnfile)
		if ids != o_ids: # sequence order has changed... revert to original order
			reorderFastaByLabel(alnfile, o_ids, ids, seqs)
			ids,seqs,need_aln = parseFasta(alnfile)
		
		assert o_ids==ids	
		assert not need_aln
	else:
		alnfile = msafilename
	
	seq2msa=[]
	for i in xrange(len(ids)):
		seq2msa.append(mapSeq2MSA(seqs[i]))
	
	return ids, seqs, seq2msa
def getSequencesAndMap(msafilename):
	# if all sequences have same length, no alignment necessary
	o_ids, o_seqs, need_aln = parseFasta(msafilename)
	if need_aln: 
		print "ERROR: sequences not aligned!"
		print len(o_ids),len(o_seqs)
		for i in xrange(len(o_ids)):
			print o_seqs[i]+"\t"+o_ids[i]
		sys.exit()
	alnfile = msafilename
	
	seq2msa=[]
	for i in xrange(len(o_ids)):
		seq2msa.append(mapSeq2MSA(o_seqs[i]))
	
	return o_ids, o_seqs, seq2msa

def remapResPos4Atoms(atoms, map2msa, msa2seq, seq):# all sequence positions start at 1
	newatoms= copy.deepcopy(atoms)
	for a in atoms:
		tmp = atoms[a]["rnum"]
		if tmp in map2msa:
			s2m = map2msa[ tmp ]
			if s2m in msa2seq:
				m2s = msa2seq[s2m]
				newatoms[a]["rnum"] = m2s
				newatoms[a]["rtype"] = AA_1_to_3[ seq[m2s-1] ] 
	return newatoms
def remapAtoms2MSA(atoms, map2msa):# all sequence positions start at 1
	newatoms= copy.deepcopy(atoms)
	for a in atoms:
		tmp = atoms[a]["rnum"]
		if tmp in map2msa:
			s2m = map2msa[ tmp ]
			newatoms[a]["rnum"] = s2m
			#newatoms[a]["rtype"] = AA_1_to_3[ seq[m2s-1] ] 
	return newatoms

def writeIdsAndSeqs2Fasta(filename, labels, sequences):
	
	assert len(labels)==len(sequences)
	f=open(filename,"w")
	for i in xrange(len(labels)):
		f.write(">%s\n%s\n"%(labels[i],sequences[i]))
	f.close()

def degap(string):
	return string.replace("-","")
def mergeSeqs2FastaAndGetIndices(seqfile, newlabels, newseqs, newfilename):
	assert len(newlabels)== len(newseqs)
	newseq_i = [ ]
	
	labels, seqs, need_aln = parseFasta(seqfile)
	assert len(seqs) == len(labels)	
	
	for i in xrange(len(newseqs)):
		n = degap(newseqs[i])
		found = None
		for si in xrange(len(seqs)):
			s = degap(seqs[si])
			if n==s:
				found = si
				break
		if found != None:
			newseq_i.append(found)
		else:
			seqs.append(degap(newseqs[i]))
			labels.append(newlabels[i])
			newseq_i.append(len(seqs)-1)
	
	
	assert len(newseq_i) == len(newseqs)
	
	writeIdsAndSeqs2Fasta(newfilename, labels, seqs)
			
	return newseq_i

def seqIdentity(s1,s2):
	l=min(len(s1),len(s2))
	m=0.0
	for i in xrange(l):
		if s1[i] == s2[i]:
			m += 1.0
	return m/l
	

def getContactInDegreePerResidue(contacts):
	indegree={}
	for c in contacts:
		res1 = c[0]
		res2 = c[1]
		dists = c[2]
		if res1 in indegree:
			indegree[res1] +=1
		else:
			indegree[res1] = 1
		
		if res2 in indegree:
			indegree[res2] +=1
		else:
			indegree[res2] = 1
	
	return indegree

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


def findClosestSeqInList(seq, seqlist):
	minHD = len(seq)
	minHD_i = -1
	for si in xrange(len(seqlist)):
		s = seqlist[si]
		h = hamDist(seq,s)
		if h<minHD: 
			minHD=h
			minHD_i = si
	return minHD_i, minHD, seqlist[minHD_i]
def hamDist(seq1, seq2):
	"""@sig void hamDist(String s1, String S2)"""
	#this elegant implementation was found on wikipedia "Hamming distance",2008-09-17
	#assert len(s1) == len(s2)
	#return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
	if len(seq1) != len(seq2):
		print "ERROR: unequal lengths"#
		return None
	count=0
	for i in xrange(len(seq1)):
		if seq1[i] != seq2[i]:
			count+=1
	return count

# MAIN PROGRAM

if __name__=="__main__":
	# SETTING UP
	
	if True: # HANDLING COMMAND LINE OPTIONS...
		##################
		# HANDLE OPTIONS (-h for help)
		description="%prog prepares native contact files in XML format for PROFASI as modified by TS.\n Default values for options are given in parentheses."
		version="2015-09-17; written by Tobias Sikosek"
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
		op.add_option('-l','--MSA',action='store',type='string',dest='MSA_file',help="filename of Multiple Sequence Alignment in a format read by BioPython. Will generate compatible contact lists for each aligned sequence.")
		op.add_option('-v','--verbosity',action="store",type="int",dest="verbosity",help="0: minimal output, 1: important output, 2: debugging")
		op.add_option('-C','--nchains',action="store",type="int",dest="nchains",help="number of identical chains")
		op.add_option('-H','--removeSharedFrom', action="store", type="string", dest="removeSharedFrom", help="PROFASI contact file. All contacts present in this file will be removed from the output contacts.")
		op.add_option('-i','--intersectWith', action="store", type="string", dest="intersectWith", help="PROFASI contact file. All contacts NOT present in this file AND the regular input file (-s) will be removed from the output contacts.")
		
		op.add_option('-D', '--homodimer', action='store', type="string", dest="dimerfile", help="PDB file containing a homodimer of the protein")
		op.add_option('-I','--indegree',action='store_true',dest="calcContactIndegree",help='calculate indegree of amino acids in contact network (i.e. how many contacts per amino acid)')
		op.add_option('-f','--chain',action="store",type="string",dest="defaultChain",help="single chain to use. defaults to A")
		
		
		op.set_defaults(inPDBFilename="",
						contactType="RESIDUE",
						contactAtoms="CA",
						chainSep=3,
						distCutoff=6.0,
						wellDepth=1.0,
						wellWidth=0.5,
						repulDist=0.0,
						repulSteep=1.0,
						ensembleCutoff=1.0,
						normEnergy=0,
						gmxpath="",
						scmpath="",
						label="",
						seqFrom=0,
						includeNonNative=False,
						selectedRange="",
						smoothGaussian=False,
						MSA_file="",
						verbosity=1,
						nchains=1,
						removeSharedFrom = "",
						intersectWith = "",
						dimerfile="",
						calcContactIndegree=False,
						defaultChain='A'
						)
	
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
		verbosity = opt.verbosity
		MSA_file = opt.MSA_file
		nchains = opt.nchains
		removeSharedFrom = opt.removeSharedFrom
		intersectWith = opt.intersectWith
		dimerfile = opt.dimerfile
		calcContactIndegree = opt.calcContactIndegree
		defaultChain = opt.defaultChain
	
		assert contactType in ["CUTOFF","SHADOW","HYBRID","RESIDUE"]
		assert distCutoff > 1.0, "cutoff in Angstrom!"
	
		restraintType = "native_restraints"
		assert restraintType in ["distance_restraints", "native_restraints"]
	
		contactPotential = "FMULTIGAUSS"
		if smoothGaussian:
			contactPotential = "FMULTIGSMOOTH"
		assert contactPotential in ["FMULTIGAUSS","FMULTIGSMOOTH"]
		if verbosity >=1: print "Contact potential is: %s "%contactPotential
		assert inPDBFilename != ""
		if normEnergy > 0 and verbosity >= 1:
			print "Normalizing contact energy per native structure to:", normEnergy
		
		removed_string = ""
		
	if True: # MISC HIDDEN OPTIONS:
		plotDimerPotentials = False
		plotDimerContactMaps = False
		
		plotFinalPotentials = False
		plotFinalContactMaps = True
		
		if any([plotDimerPotentials,plotDimerContactMaps,plotFinalPotentials,plotFinalContactMaps]):
			import matplotlib.pyplot as p
			import matplotlib.cm as cm
		
		figdpi = 300
	
	if True: # SETTING UP RESIDUE RANGES... IF APPLICABLE
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
	
		
	
	
	
	
	if True: # WHICH ATOM TYPES ARE USED FOR CONTACTS ?
		contactAtoms = contactAtoms.split(",")
	
		if contactType == "SHADOW" and contactAtoms != ['CA']:
			if verbosity >= 1: print "SHADOW only allows CA atoms!"
			contactAtoms = ['CA']
	
		if label=="":
			label = "_"+inPDBFilename[:-4]
	
		if verbosity >= 1 and contactType=="CUTOFF": print "Atoms to consider for CUTOFF contacts:",contactAtoms
	
		assert all([ a in BB_atoms for a in contactAtoms])
	
		if contactType in ["HYBRID","SHADOW"]:
			if verbosity >= 1: print "Also considering CA atoms for SHADOW contacts."
	
	
	
	
	# OBTAIN CONTACTS
	
	if True: # obtain multiple structure models from single PDB file (-s option)
		# setting up lists for PDB models
		atoms_models = []
		models_contact_lists = []
		models_scenergy_lists = []
		model_sequences = []
		model_labels = []
		model_chains = []
		
		models = splitModels(inPDBFilename, model_labels)
		if verbosity >=1: print "\n\n%i model(s) found in %s."%(len(models), inPDBFilename)
		
	
	
	
	for m in models: # Iterate over all models (from -s option), collect contacts for each
		contacts = []
		con_c=[]
		con_s=[]
		
		atoms = atomsFromPDB(m)
		seq = getSequence(atoms)
		model_sequences.append( seq )
		chainlabels = getChainIDs(atoms)
		if verbosity >=2: print "model %s chains:"%m, chainlabels
		model_chains.append(chainlabels)
		
		if verbosity >=2: print "model %s sequence:"%m, seq
		atoms_models.append(atoms)
		
		if verbosity >=2: print "num atoms",len(atoms), m
		
		if contactType=="CUTOFF":
			m_contacts = getCutoffContacts(atoms, atomtypes=contactAtoms, cutoff=distCutoff, chainDist=chainSep, verbose=False, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		elif contactType=="SHADOW":
			m_contacts = getShadowContacts(pdbfilename=m, outfilename=inPDBFilename[:-4]+".shadow", cutoff=distCutoff, chainDist=3, gmxpath=gmxpath, scmpath=scmpath, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		elif contactType=="HYBRID":
			con_c = getCutoffContacts(atoms, atomtypes=contactAtoms, cutoff=distCutoff, chainDist=chainSep, verbose=False, res2SS = None, selection=selection, selectExclusive=selectExclusive)
			con_s = getShadowContacts(pdbfilename=m, outfilename=inPDBFilename[:-4]+".shadow", cutoff=distCutoff, chainDist=3, gmxpath=gmxpath, scmpath=scmpath, res2SS = None, selection=selection, selectExclusive=selectExclusive)
			m_contacts = mixCutoffShadow(con_c, con_s)
		elif contactType=="RESIDUE":
			if len(chainlabels)>1 and dimerfile == "":
				m_contacts = getInterchainResidueContacts(atoms, chainSep, distCutoff, verbose=False, selection=selection, selectExclusive=selectExclusive, interChainContacts=True, intraChainContacts=True)
			else:
				m_contacts = getInterchainResidueContacts(atoms, chainSep, distCutoff, verbose=False, selection=selection, selectExclusive=selectExclusive, interChainContacts=False, intraChainContacts=True)#getResidueContacts(atoms, cutoff=distCutoff, chainDist=chainSep,verbose=False, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		else:
			print "unknown contact type: ",contactType
			sys.exit(1)
		
		models_contact_lists.append(m_contacts)
	
	
	
	
	if MSA_file != "": # use MSA (option -l) for filtering consensus contacts; residue indices now refer to column indices in MSA !
		
		if False:# include PDB models (-i option) into alignment?
			# include pdb sequences in MSA
			newseqfilename = os.path.splitext(MSA_file)[0] + "_%s"%(os.path.splitext(inPDBFilename)[0]) + os.path.splitext(MSA_file)[1]
			if len(model_labels) != len(model_sequences):
				model_labels = [inPDBFilename[:-4]+"_m%s"%(ms+1) for ms in xrange(len(model_sequences))]
				print model_labels
		
			#mergeSeqs2FastaAndGetIndices(MSA_file, model_labels, model_sequences, newseqfilename)
		
			#if verbosity >=1: print "Wrote new fasta file (not aligned yet):",newseqfilename
		
			# re-align all sequences
			msa_labels, msa_seqs, seq2msa = alignSequencesAndMap(newseqfilename)
			
			if verbosity >=1: print "Aligned sequences written to fasta file: ",newseqfilename
		else:
			msa_labels, msa_seqs, seq2msa = getSequencesAndMap(MSA_file)
			model_labels=[]
		print msa_labels, len(msa_labels)
		modelSeqIds = [None for i in xrange(len(model_sequences))]
		
		if verbosity >=1: print "\nTHE ALIGNMENT:"
		if verbosity >=1: print sequencePosMarker(len(msa_seqs[0])) 
		
		# the PDB model sequences from alignemnt,i.e. with gaps
		model_alignedSeqs=[None for i in xrange(len(model_sequences))] 
		
		# print sequences from alignment and collect PDB sequences along the way (marked with "*")
		for i in xrange(len(msa_labels)):
			models_found_counter = 0
			for j in xrange(len(model_sequences)):
				if degap(msa_seqs[i]) == model_sequences[j] and "pdb" in msa_labels[i]:
					modelSeqIds[j] = i
					model_alignedSeqs[j] = msa_seqs[i]
					model_labels.append(msa_labels[i])
					models_found_counter += 1
			if verbosity >=1: 
				if models_found_counter ==0:
					print ">%s\n%s"%(msa_labels[i],msa_seqs[i])
				else:
					print ">%s (%i models)\n%s"%(msa_labels[i],models_found_counter,msa_seqs[i])
		
		if verbosity >=1: print sequencePosMarker(len(msa_seqs[0])) 
		
		#print modelSeqIds, model_alignedSeqs
		
		for i in xrange(len(model_alignedSeqs)):
			if model_alignedSeqs[i] == None:
				print "ERROR: could not find sequence from model in alignment:",model_sequences[i]
				#print "closest sequence:"
				#print findClosestSeqInList(model_sequences[i], [degap(mi) for mi in msa_seqs])
		
		assert not None in model_alignedSeqs
		
		# obtain consensus contacts (from all PDB models from option -s) with residue numbering from MSA
		#print modelSeqIds
		finalcontacts = filterSuperEnsembleContacts(models_contact_lists, 
													atoms_models, 
													model_chains, 
													contactAtoms, 
													ensembleCutoff, 
													contactType, 
													map2msa=seq2msa, 
													modelSeqMsaIDs=modelSeqIds
													)
	
	else: # filter consensus contacts, assuming all input models have the same sequence length
		# obtain consensus contacts (from all PDB models from option -s) - assuming all sequences have the same length !
		finalcontacts = filterSuperEnsembleContacts(models_contact_lists, 
													atoms_models, 
													model_chains, 
													contactAtoms, 
													ensembleCutoff, 
													contactType
													)
	
	
	
	
	
	
	if dimerfile != "": # get dimer contacts from file (option -D)
		dimer_model_labels = []
		generic_dimer_model_labels = [] #deprecated: automatically generated labels based on filename, not very informative. MSA labels should be preferred (see dimer_model_labels)
		dimer_models_contact_lists = []
		dimer_model_sequences = []
		dimer_atoms_models = []
		dimer_model_chains = []
		dimer_models = splitModels(dimerfile, generic_dimer_model_labels)
		
		
		if verbosity >=1: print "\n\n DIMERS: %i model(s) found in %s."%(len(dimer_models), dimerfile)
		
		for m in dimer_models:
			atoms = atomsFromPDB(m)
			dimer_atoms_models.append(atoms)
			seq = getChainSequences(atoms)
			dimer_model_sequences.append( seq )
			dimer_chainlabels = getChainIDs(atoms)
			if verbosity >=2: print "dimer model %s chains:"%m, dimer_chainlabels
			dimer_model_chains.append(dimer_chainlabels)
			
			if contactType=="RESIDUE":
				dm_contacts = getInterchainResidueContacts(atoms, chainSep, distCutoff, verbose=False, selection=selection, selectExclusive=selectExclusive, interChainContacts=True, intraChainContacts = False)
			else:
				print "ONLY 'RESIDUE' contacts supported for dimers!"
				sys.exit(1)
		
			#print dm_contacts, len(dm_contacts)
			dimer_models_contact_lists.append(dm_contacts)
		
		
		if MSA_file != "":
			dimer_model_msa_indices  =[]
			dimer_model_representative_seqs = []
			for a in xrange(len(dimer_models)):
				dimer_aln_label = None
				dimer_aln_index = None
				
				min_dms = -1
				minlength = 10000
				for di in xrange(len(dimer_model_sequences[a])):
					d = dimer_model_sequences[a][di]
					if len(d)<minlength:
						minlength = len(d)
						min_dms = di
				#print a,len(dimer_model_sequences)
				dimer_sequence = dimer_model_sequences[a][min_dms]
				if verbosity >= 2: print "Selected dimer chain sequence:",dimer_sequence, " - Chain:", dimer_model_chains[a][min_dms]
				for mi in xrange(len(msa_seqs)):
					if dimer_sequence == degap(msa_seqs[mi]):
						dimer_aln_index = mi
						dimer_aln_label = msa_labels[mi]
						#print "found in MSA:",dimer_aln_label, dimer_aln_index
						break
				
				#print "Dimer chain sequences:",dimer_model_sequences[a], " - Selected:", max_dms
			
				if dimer_aln_label != None:
					if verbosity >= 2: print "Found dimer sequence in alignment with label:", dimer_aln_label
					dimer_model_labels.append(dimer_aln_label)
					dimer_model_representative_seqs.append(dimer_sequence)
					dimer_model_msa_indices.append(dimer_aln_index)
				else:
					print "ERROR: could not find dimer sequence in alignment!"
					print a, dimer_model_sequences[a]
					sys.exit(1)
			
			dimer_finalcontacts = filterSuperEnsembleContacts(dimer_models_contact_lists, 
																dimer_atoms_models, 
																dimer_model_chains, 
																contactAtoms, 
																ensembleCutoff, 
																contactType, 
																map2msa=seq2msa, 
																modelSeqMsaIDs=dimer_model_msa_indices,
																isDimer=True
																)
			
			
			if len(dimer_model_labels) == 0:
				dimer_model_labels = generic_dimer_model_labels
			
			assert len(dimer_model_labels)==len(dimer_models_contact_lists)==len(dimer_model_sequences)==len(dimer_atoms_models)==len(dimer_model_chains)==len(dimer_models), map(str,[len(dimer_model_labels),len(dimer_models_contact_lists),len(dimer_model_sequences),len(dimer_atoms_models),len(dimer_model_chains),len(dimer_models)])
			
			print "Dimer summary:"
			for i in xrange(len(dimer_models)):
				print i, dimer_model_labels[i], str(len(dimer_atoms_models[i]))+ " atoms ",str(len(dimer_models_contact_lists[i]))+"("+str(len(dimer_finalcontacts))+") contacts(final) ", len(dimer_model_sequences[i]),"seqs ", "chains:"+str(dimer_model_chains[i]), "shorter seq:",dimer_model_representative_seqs[i]
			
			
			
		else:
			dimer_model_labels = [] #TODO: complete this for the case that dimers are used without MSA?
			dimer_finalcontacts = filterSuperEnsembleContacts(dimer_models_contact_lists, 
																dimer_atoms_models, 
																dimer_model_chains, 
																contactAtoms, 
																ensembleCutoff, 
																contactType
																)
	
	
	
	
	
	
	if removeSharedFrom != "":# filter contacts from XML input file (option: -H)
		rem = readXMLContacts(removeSharedFrom)
		print rem.keys(), len(rem)
		toBeRemoved=[]
		for c in xrange(len(finalcontacts)):
			r1 = finalcontacts[c][0]
			r2 = finalcontacts[c][1]
			print "removing contacts from file %s:"%removeSharedFrom
			if (r1,r2) in rem:
				print (r1,r2)
				toBeRemoved.append(c)
		print len(toBeRemoved), "contacts removed"
		removed_string = "_%iremShared"%len(toBeRemoved)
		for tbr in sorted(toBeRemoved, reverse=True):
			del finalcontacts[tbr]
	
	if intersectWith != "":
		rem = readXMLContacts(intersectWith)
		print rem.keys(), len(rem)
		toBeRemoved=[]
		for c in xrange(len(finalcontacts)):
			r1 = finalcontacts[c][0]
			r2 = finalcontacts[c][1]
			print "removing contacts from file %s:"%removeSharedFrom
			if not (r1,r2) in rem:
				print (r1,r2)
				toBeRemoved.append(c)
		print len(toBeRemoved), "contacts removed"
		removed_string = "_%iremUnique"%len(toBeRemoved)
		for tbr in sorted(toBeRemoved, reverse=True):
			del finalcontacts[tbr]
	
	if True:# generate automatic output file name
		consStatement = "_cons%s"%str(ensembleCutoff)
		outatoms = "".join([i for i in contactAtoms])
		width_label = "_w"+str(round(wellWidth,1))
		repul_label = "_R"+str(round(repulDist,1))
		nonnat_label = "_x"+str(includeNonNative)
		outfilename = "%s_%s_%s_c%s_n%s%s%s_%s%s%s_%s%s"%(inPDBFilename.split("/")[-1],contactPotential,contactType,str(distCutoff),str(chainSep),width_label,repul_label,outatoms,consStatement,nonnat_label,selectedRange,removed_string)
	
	
	
	
	
	if contactType=="HYBRID" and not 'CA' in contactAtoms:
		contactAtoms.append('CA')
	
	
	
	pbond_legend = None
	
	if MSA_file != "": # write contacts to file for each sequence in the MSA
		#msa_outfilenames = []
		
		#output filename also contains the name of the sequence from the MSA
		msa_outfilename = "%s_%s_%s_c%s_n%s%s%s_%s%s%s_%s_%s_nch%i"%( inPDBFilename.split("/")[-1], contactPotential, contactType, str(distCutoff), str(chainSep), width_label, repul_label, outatoms, consStatement, nonnat_label, selectedRange, "MSA-POSITIONS",nchains )
		# Writing PROFASI XML file:
		# POSITIONS CORRESPOND TO MSA !!! 
		msa_outfilename, msa_co = writeProfasiContactsFMULTIGAUSS(finalcontacts, 
															atoms, 
															contactAtoms, 
															msa_outfilename, 
															repulDist, 
															repulSteep, 
															width=wellWidth, 
															depth=wellDepth, 
															norm=normEnergy, 
															label=label, 
															includeNonNative=includeNonNative, 
															chainSep=chainSep, 
															smoothGaussian=smoothGaussian, 
															restraintType=restraintType, 
															nchains=nchains, 
															genericSequenceOfLength=len(msa_seqs[0]) # number of columns in MSA
															)
		if verbosity >=2: print "wrote Profasi XML contacts (%s):"%mlab,msa_outfilename, " - CO =",msa_co
		
		
		
		for mi in xrange(len(msa_seqs)): # for each sequence in the alignment ...
			mseq = msa_seqs[mi]
			#print mseq
			mlab = msa_labels[mi]
			msa2seq = mapMSA2Seq(mseq) # map positions from aligned sequence to the same sequence without gaps
			
			#which PDB sequence (from -s input) most closely matches this aligned sequence?
			closestModel2Seq = -1
			dist = -1.0
			for modi in xrange(len(model_alignedSeqs)):
				tmp = seqIdentity(mseq,model_alignedSeqs[modi])
				if tmp >dist:
					dist = tmp
					closestModel2Seq = modi
			
			if verbosity >=2: 
				print "Closest PDB model to sequence named '%s' has index %i and is named '%s'. Its sequence is:\n%s"%( mlab, closestModel2Seq, model_labels[closestModel2Seq], model_alignedSeqs[closestModel2Seq] )
				print "Indices of PDB model sequences in the alignment: ",modelSeqIds
			
			#map_model2mseq = mapMSASeq1Seq2(model_alignedSeqs[closestModel2Seq], mseq)
			# this ONLY makes sense for chimera pseudobonds, where the displayed PDB protein is 'closestModel2Seq'
			map_model2mseq_chimera = mapMSASeq1ToDegappedSeq2(model_alignedSeqs[closestModel2Seq], mseq)
			
			# Change residue numbers and names of atom coordinates in PDB model. 
			# From ungapped PDB sequence to MSA positions
			msa_atoms = remapAtoms2MSA( atoms_models[closestModel2Seq], seq2msa[ modelSeqIds[closestModel2Seq] ])
			
			if verbosity >=2: print "Mapping contacts to sequence " + mlab
			if verbosity >=2: print mseq + "\n" + sequencePosMarker(len(mseq))
			
			# different contact files are written for each sequence in the MSA
			
			#output filename also contains the name of the sequence from the MSA
			msa_outfilename = "%s_%s_%s_c%s_n%s%s%s_%s%s%s_%s_%s_nch%i"%( inPDBFilename.split("/")[-1], contactPotential, contactType, str(distCutoff), str(chainSep), width_label, repul_label, outatoms, consStatement, nonnat_label, selectedRange, mlab.replace("/",""),nchains )
			
			# Writing PROFASI XML file:
			# Positions should match MSA sequence, will be converted to ungapped sequence position just before writing XML file
			msa_outfilename, msa_co = writeProfasiContactsFMULTIGAUSS(finalcontacts, 
																msa_atoms, 
																contactAtoms, 
																msa_outfilename, 
																repulDist, 
																repulSteep, 
																width=wellWidth, 
																depth=wellDepth, 
																norm=normEnergy, 
																label=label, 
																includeNonNative=includeNonNative, 
																chainSep=chainSep, 
																smoothGaussian=smoothGaussian, 
																restraintType=restraintType, 
																mapMSA2Seq=msa2seq, 
																nchains=nchains, 
																generateChainInfo = True, 
																aligned_sequence=mseq
																)
			if verbosity >=2: print "wrote Profasi XML contacts (%s):"%mlab,msa_outfilename, " - CO =",msa_co
			
			if mlab in model_labels:
				# Writing Chimera pseudobond file:
				msa_pbond_filename = msa_outfilename[:-4]+".pseudobonds"
				writeChimeraPseudobonds(finalcontacts, 
										msa_pbond_filename, 
										legend=pbond_legend, 
										mapMSA2Seq=map_model2mseq_chimera, 
										verifyPDBatoms=msa_atoms,
										nchains=nchains, 
										generateChainInfo = True
										) #atoms_models[closestModel2Seq])
				if verbosity >=2: print "wrote Chimera pseudobonds (MSA):",msa_pbond_filename
			
			if dimerfile != "":
				if mlab in dimer_model_labels:
					#print "corresponds to dimer:",mlab,,
					dimer_pbond_filename = msa_outfilename[:-4]+"_DIMER.pseudobonds"
					writeChimeraPseudobonds(dimer_finalcontacts, 
											dimer_pbond_filename, 
											legend=pbond_legend, 
											mapMSA2Seq=map_model2mseq_chimera, 
											verifyPDBatoms=msa_atoms
											)
					if verbosity >=2: print "wrote Chimera pseudobonds (dimer):",dimer_pbond_filename
				
				#print dimer_finalcontacts
				dimer_msa_outfilename, d_co = writeProfasiContactsFMULTIGAUSS(dimer_finalcontacts, 
																		msa_atoms, 
																		contactAtoms, 
																		msa_outfilename+"_DIMER", 
																		repulDist, 
																		repulSteep, 
																		width=wellWidth, 
																		depth=wellDepth, 
																		norm=normEnergy, 
																		label=label+"_DIMER", 
																		includeNonNative=includeNonNative, 
																		chainSep=chainSep, 
																		smoothGaussian=smoothGaussian, 
																		restraintType=restraintType, 
																		mapMSA2Seq=msa2seq, 
																		nchains=nchains, 
																		aligned_sequence=mseq
																		)
				if plotDimerPotentials:
					plotMultiGaussianContacts(dimer_finalcontacts,1.0,0.5,dimer_msa_outfilename[:-4]+"_potentials.png",smoothGaussian=smoothGaussian,dpi=figdpi)
				
				if plotDimerContactMaps:
					plotContactMap(dimer_finalcontacts, mlab, distCutoff, chainSep, msa_atoms, contactType, fname=dimer_msa_outfilename[:-4]+"_contacts.png",dpi=figdpi)
				
				if verbosity >=2: print "wrote Profasi XML contacts (dimer):",dimer_msa_outfilename, " - CO =",d_co
	else: # write contacts to file
		if verbosity >=1: print "Sequence from model",seqFrom
		# Writing PROFASI XML file:
		outfilename, co = writeProfasiContactsFMULTIGAUSS(finalcontacts, 
														atoms_models[seqFrom], 
														contactAtoms, 
														outfilename, 
														repulDist, 
														repulSteep, 
														width=wellWidth, 
														depth=wellDepth,
														norm=normEnergy,
														label=label,
														includeNonNative=includeNonNative,
														chainSep=chainSep,
														smoothGaussian=smoothGaussian,
														restraintType=restraintType,
														nchains=nchains
														)
		if verbosity >=2: print "wrote Profasi XML contacts:",outfilename, " - CO =",co
		
		# Writing Chimera pseudobond file:
		writeChimeraPseudobonds(finalcontacts, 
								outfilename[:-4]+".pseudobonds",
								legend=pbond_legend
								)
		if verbosity >=2: print "wrote Chimera pseudobonds:",outfilename[:-4]+".pseudobonds"
		if verbosity >= 1: print "Total native energy:",-len(finalcontacts)*wellDepth
	
	
	
	
	if verbosity >= 1: print "Total native energy:",-len(finalcontacts)*wellDepth
	if calcContactIndegree: # calculate contact in-degree for each residue
		residues = getResidueDict(atoms_models[seqFrom])
		indegfile = open( outfilename[:-4]+".ncontacts","w")
		indeg = getContactInDegreePerResidue(finalcontacts)
		print model_sequences[seqFrom]
		print residues.keys()
		nskip=0
		for k in xrange(len(model_sequences[seqFrom])+1):
			if k+1 in indeg and k+1 in residues:
				deg = indeg[k+1]
			elif not k+1 in indeg and k+1 in residues:
				deg = 0
			else:
				nskip +=1
				continue
			if verbosity >= 1: print k+1, model_sequences[seqFrom][k-nskip], deg
			indegfile.write(",".join(map(str,[k+1, model_sequences[seqFrom][k-nskip], deg]))+"\n")
		indegfile.close()
	
	
	
	if plotFinalPotentials:# plot graph showing potential functions for all contacts
		plotMultiGaussianContacts(finalcontacts,1.0,0.5,outfilename[:-4]+"_potentials.png",smoothGaussian=smoothGaussian,dpi=figdpi)
	
	if plotFinalContactMaps:# plot consensus contact map
		if MSA_file != "":
			useMSAPos = len(msa_seqs[0])
		else:
			useMSAPos = None
		plotContactMap(finalcontacts, inPDBFilename, distCutoff, chainSep, atoms, contactType, fname=outfilename[:-4]+"_contacts.png", dpi=figdpi, useMSAPos=useMSAPos)
	
