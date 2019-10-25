import sys,math,os,random,subprocess,re,numpy


def parseDataFromFolders(nrep, temps, minsteps, maxsteps):
	maxtime = 0
	
	
	etotPerTemp = [[] for t in xrange(len(temps))]
	
	for r in xrange(nrep):
		for line in open( "n%i/rt"%r).readlines():
			tmp = line.replace("\x00","").strip().split()
		
			
			
	
			#if len(tmp) != nDataCols:
			#	continue
			#print tmp
			if tmp[0]=='': print tmp
			try:
				time = int(tmp[0])
			except ValueError:
				print tmp
				continue
		
			if minsteps > 0 and time < minsteps:
				continue
			if maxsteps > 0 and time > maxsteps:
				continue
		
			if time > maxtime: maxtime = time
		
		
			Ti = int(tmp[1])
			
				
		
			Etot= float(tmp[2])
		
		
			etotPerTemp[Ti].append(Etot)
		
	if maxtime < minsteps:
		print "ERROR: no data - lower minsteps!"
		sys.exit(1)
	
	return etotPerTemp
def getTemps(filename):
	t = open(filename)
	h = t.readline().strip("#").strip().split()
	d = {}
	for line in t.readlines():
		tmp = line.strip().split()
		index = int(tmp[0])
		Tprf = float(tmp[1])
		Tkel = float(tmp[2])
		beta = float(tmp[3])
		d[index] = (Tprf,Tkel,beta)
	t.close()
	return d

def loadModulesScinet():
	cmd = "module load gcc/4.8.1  openmpi/gcc/1.6.4"
	s = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	#print cmd
	out = s.stdout.read()
	err = s.stderr.read()
	#print out
	if err.strip() != "": print err
def extractPDB(r, time, label, outputdir, PRFLOC="~/PROFASI/app/bin/"):
	fname = "%s_n%i_t%i.pdb"%(label, r, time)
	if not os.path.exists("%s/%s"%(outputdir,fname)):
		cmd = "module load gcc/4.8.1  openmpi/gcc/1.6.4; %sextract_snapshot -o %s/%s n%i/traj -c %i -f pdb"%(PRFLOC, outputdir, fname, r, time)
		#cmd = "%sextract_snapshot -o %s/%s n%i/traj -c %i -f pdb"%(PRFLOC, outputdir, fname, r, time)
		s = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		#print cmd
		out = s.stdout.read()
		err = s.stderr.read()
		#print out
		if err.strip() != "": print err
	
	
	return fname

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
						contacts.append( (i,j,natoms,mindist) )
					else:
						contacts.append( (i,j) )
					#print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	
	
	return contacts
def getCv(temps, etotPerTemp, maxTmCut=2000, polyfitDegree=None):
	#k_B=0.0019872041
	#kBprf = k_B / kcal_per_mol_in_prf_energy * kelvin_in_pru
	tkel = [temps[t][1] for t in reversed(xrange(len(temps)))]
	tprf = [temps[t][0] for t in reversed(xrange(len(temps)))]
	cvs = [ (1.0/math.pow(tprf[t],2)) * ( numpy.mean( [math.pow(i,2) for i in etotPerTemp[t]] ) - math.pow(numpy.mean([i for i in etotPerTemp[t]]),2)    )   for t in reversed(xrange(len(temps))) ]
	tmp = [cvs[c] for c in xrange(len(cvs)) if tkel[c] <= maxTmCut ]
	
	
	if polyfitDegree != None:
		
		fit = numpy.polyfit(tkel,cvs,polyfitDegree)
		#print fit
		finegrid = numpy.arange(min(tkel),max(tkel),(max(tkel)-min(tkel)) / 100.0) 
		#print finegrid
		fineCv = numpy.polyval(fit, finegrid)
		
		maxcv = max(fineCv)
		maxcv_i = list(fineCv).index(maxcv)
		Tmelti = maxcv_i
		Tmelt = finegrid[Tmelti]
		#print "T grid:",finegrid, Tmelti,Tmelt
		#print "Cv grid:", fineCv, maxcv, maxcv_i
		cvs = list(fineCv)
	else:
		Tmelt = None
		Tmelti = None
		maxcv  = max(tmp)
		for i in xrange(len(cvs)):
			if cvs[i]==maxcv:
				Tmelt = tkel[i]
				Tmelti = i
	
	
	return cvs, Tmelt, Tmelti, maxcv
def atomDist(a1,a2,atoms):
	d = math.sqrt( math.pow(atoms[a1]["x"] - atoms[a2]["x"], 2)
					+ math.pow(atoms[a1]["y"] - atoms[a2]["y"], 2)
					+ math.pow(atoms[a1]["z"] - atoms[a2]["z"], 2)
					)
	#print d
	return d

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
	for line in open(filename).readlines():
		if line[:4]=="ATOM":
			a = parseATOMline(line)
			atoms[a["anum"]] = a
	return atoms







def plotContactMap(n, contacts, contactsA, contactsB, inPDBFilename, cutoff, chainSep, atoms, contactType, figsavename):
	
	
	minP = min([contacts[i] for i in contacts.keys()])
	maxP = max([contacts[i] for i in contacts.keys()])
	
	matrix = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	#matrix1 = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	#matrix2 = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	
	for i in xrange(n):
		for j in xrange(n):
			if (i,j) in contacts:
				matrix[i][j] = contacts[(i,j)]
				#print i,j,contacts[(i,j)]
			if (i,j) in contactsA and not (i,j) in contactsB:
				matrix[j][i] = maxP*0.1
			elif (i,j) in contactsB and not (i,j) in contactsA:
				matrix[j][i] = maxP*1.0
			elif (i,j) in contactsA and (i,j) in contactsB:
				matrix[j][i] = maxP*0.4
	
	title = "%s\ncutoff=%s;csep=%i;n=%i;%s"%(inPDBFilename,str(round(cutoff,2)),chainSep,len(contacts),contactType)
	
	
	fig = p.figure(figsize=(5,5))
	ax = p.subplot(111)
	
	imax = p.matshow(matrix, cmap=cm.binary, aspect='equal')
	
	majorLocator   = MultipleLocator(10)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator   = MultipleLocator(1)
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(majorFormatter)
	ax.xaxis.set_minor_locator(minorLocator)
	ax.yaxis.set_major_locator(majorLocator)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.yaxis.set_minor_locator(minorLocator)
	#fig, ax = plt.subplots()
	p.title(title)
	cbar = p.colorbar()
	p.xlabel("residue position")
	p.ylabel("residue position")
	cbar.set_label("P(contact)")
	#imax.get_axes().set_axisbelow(True)
	#p.tight_layout()
	p.savefig(figsavename)
	#p.show()

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
			contacts[(r1,r2)] = tmp
	return contacts


def writeContactHisto(n, histo, filename):
	f = open(filename,"w")
	f.write(str(n) + "\n")
	f.write(str(histo) + "\n")
	f.close()
	return True


def getRTkeys():
	return [i.strip() for i in open("n0/rtkey").readlines()[1:]]
	
if __name__=="__main__":
	# module load intel/15.0.1 python/2.7.8
	
	print "available temperatures:\nindex (prf,kelvin,beta)"
	temps = getTemps("n0/temperature.info")
	for i in xrange(len(temps)): print i,temps[i]
	
	mode = sys.argv[1]
	assert mode in ["foldedA","foldedB","unfolded","transA","transB","all"]
	#mode = "foldedA"
	#mode = "foldedB"
	#mode = "unfolded"
	
	samplesize = int(sys.argv[2])
	
	EAmax = 48.8
	EBmax = 51
	print "EA=%f;EB=%f"%(EAmax,EBmax)
	
	mintime = 3000000
	maxtime = 10000000
	
	n = 32
	
	etotPerTemp = parseDataFromFolders(n, temps, mintime, maxtime)
	cvs, Tmelt, Tmelti, maxcv = getCv(temps, etotPerTemp, maxTmCut=2000, polyfitDegree=None)
	#T = int(raw_input("Select temperature:"))
	#print "chosen temperature:", T
	
	T = Tmelt
	
	print "min,max time:",mintime,maxtime
	
	ia = 11
	ib = 12
	
	print "%i replicas; data columns: %i,%i"%(n,ia,ib)
	
	cut_foldedA = 0.7
	cut_foldedB = 0.7
	cut_unfoldedA = 0.6
	cut_unfoldedB = 0.3
	
	length = 56
	
	plots = False
	
	cutoff = 6.0
	chainDist = 3
	
	
	outputdir = mode
	
	filterresidues = [45]
	
	quitAfterPDB2DATA = False
	
	#loadModulesScinet()
	
	nativeContactFile1 = "superGA.pdb_FMULTIGAUSS_RESIDUE_c6.0_n3_w0.5_R0.0_CA_cons1.0_xFalse_95contacts_E48.8.xml"
	nativeContactFile2 = "superGB.pdb_FMULTIGAUSS_RESIDUE_c6.0_n3_w0.5_R0.0_CA_cons1.0_xFalse_137contacts_E51.0.xml"
	
	natpdb1 = "2LHC_min_rmsd.pdb"
	natpdb2 = "2LHD_min_rmsd.pdb"
	
	n1_atoms = atomsFromPDB(natpdb1)
	n2_atoms = atomsFromPDB(natpdb2)
	
	native1 = readXMLContacts(nativeContactFile1)
	native2 = readXMLContacts(nativeContactFile2)
	native = native1.copy()
	native.update(native2)
	
	print len(native1.keys()),len(native2.keys()),len(native.keys()), "native contacts (native1,native2,both)"
	
	
	print "Applying same contact definition to native pdbs as simulated pdbs..."
	
	natcon1 = getResidueContacts(n1_atoms,cutoff,chainDist,includeNatoms=True)
	print "Contacts in NATIVE 1:\nres1 res2 nAtoms minDist"
	for c in natcon1: print c[0]+1, c[1]+1,c[2],c[3]
	natcon2 = getResidueContacts(n2_atoms,cutoff,chainDist,includeNatoms=True)
	print
	print "Contacts in NATIVE 2:\nres1 res2 nAtoms minDist"
	for c in natcon2: print c[0]+1, c[1]+1,c[2],c[3]
	
	label = os.path.basename(os.getcwd())
	sample = []
	samplenames = []
	
	pdb2data={}
	
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	
	datakeys = getRTkeys()
	print datakeys
	
	
	if True:
		print "\nreading data..."
		#dataPerReplica = [[] for i in n]
		data = []
		for i in xrange(n):
			sys.stdout.write("\r%s%%" % str(round(float(100.0*i)/n,2)))
			sys.stdout.flush()
			for line in open("n%i/rt"%i).readlines():
				tmp = line.strip().split()
				if len(tmp) < 4:
					continue
				#print tmp
				time = int(tmp[0].replace("\x00",""))
				if time < mintime:
					continue
				if time != -1 and time > maxtime:
					continue
				
				Ti = int(tmp[1])
				Etot = float(tmp[2])
				A = -float(tmp[ia])/EAmax
				B = -float(tmp[ib])/EBmax
				
				if ( (mode=="foldedA"  and T==Ti and A >= cut_foldedA) or 
					(mode=="foldedB"  and T==Ti and B >= cut_foldedB) or
					(mode=="transA" and T==Ti and A < cut_foldedA and A > cut_unfoldedA) or
					(mode=="transB" and T==Ti and B < cut_foldedB and B > cut_unfoldedB) or
					(mode=="unfolded" and T==Ti and A <= cut_unfoldedA and B <= cut_unfoldedB) or
					mode=="all"):
					#print A,B
					data.append( (i,tmp) )
		print
		print "%i structures available in total in state '%s'"%(len(data),outputdir)
		
		print "\nextracting pdb files..."
		#for i in xrange(samplesize):
		i=0
		while len(data)>0 and len(samplenames) < samplesize:
			sys.stdout.write("\r%s%%" % str(round(float(100.0*i)/samplesize,2)))
			sys.stdout.flush()
			i+=1
			#if len(data) == 0: 
			#	#print "no data"
			#	continue
			choice = random.choice(data)
			sample.append(choice)
			fname = extractPDB(choice[0], int(choice[1][0]), label, outputdir, PRFLOC="~/PROFASI/app/bin/")
			
			if not os.path.exists(outputdir+"/"+fname):
				print "ERROR: this file should exist, but does not:",fname
			else:
				samplenames.append(fname)
				pdb2data[fname] = choice[1]
				#print fname
			data.remove(choice)
	print
	print len(samplenames), "structures sampled of ", samplesize
	
	
	
	
	print 
	print "Writing pdb-to-data maps:"
	pdb2data_file = open(label + "_ALL_%s_N%i_pdb2data.dat"%(outputdir,samplesize),"w")
	pdblist_file = open(label + "_ALL_%s_N%i_pdblist.dat"%(outputdir,samplesize),"w")
	
	
	pdb2data_file.write(",".join(["filename","replica"]+datakeys)+"\n")
	
	for i in sorted(pdb2data.keys()):
		pdb2data_file.write(",".join([i,i.split("_")[-2][1:]]+pdb2data[i])+"\n")
		pdblist_file.write(i+"\n")
	pdb2data_file.close()
	pdblist_file.close()
	
	if quitAfterPDB2DATA:
		sys.exit()
	
	contactHistoAll = {}
	contactHistoNat1 = {}
	contactHistoNat2 = {}
	contactHistoNonNat = {}
	contactHistoFiltered = {}
	
	contactHistoAllCount = 0
	contactHistoNat1Count = 0
	contactHistoNat2Count = 0
	contactHistoNonNatCount = 0
	contactHistoFilteredCount = 0
	
	print "\nwriting histograms..."
	for i in xrange(len(samplenames)):
		sys.stdout.write("\r%s%%" % str(round(float(100.0*i)/len(samplenames),2)))
		sys.stdout.flush()
		atoms = atomsFromPDB(outputdir + "/" +samplenames[i])
		contacts = getResidueContacts(atoms, cutoff, chainDist, verbose=False, res2SS = None)
		for c in contacts:
			
			contactHistoAllCount += 1
			if c in contactHistoAll:
				contactHistoAll[c] += 1
			else:
				contactHistoAll[c] = 1
			
			if c in native1:
				contactHistoNat1Count += 1
				if c in contactHistoNat1:
					contactHistoNat1[c] += 1
				else:
					contactHistoNat1[c] = 1
			elif c in native2:
				contactHistoNat2Count += 1
				if c in contactHistoNat2:
					contactHistoNat2[c] += 1
				else:
					contactHistoNat2[c] = 1
			else:
				contactHistoNonNatCount += 1
				if c in contactHistoNonNat:
					contactHistoNonNat[c] += 1
				else:
					contactHistoNonNat[c] = 1
			
			if (c[0]+1 in filterresidues or c[1]+1 in filterresidues):
				contactHistoFilteredCount += 1
				if c in contactHistoFiltered:
					contactHistoFiltered[c] += 1
				else:
					contactHistoFiltered[c] = 1
				
	norm_contactHistoAll = {}
	norm_contactHistoNat1 = {}
	norm_contactHistoNat2 = {}
	norm_contactHistoNonNat = {}
	norm_contactHistoFiltered = {}
	
	
	
	for i in xrange(length):
		for j in xrange(length):
			if (i,j) in contactHistoAll:
				norm_contactHistoAll[(i,j)] = contactHistoAll[(i,j)] / float(contactHistoAllCount)
			if (i,j) in contactHistoNat1:
				norm_contactHistoNat1[(i,j)] = contactHistoNat1[(i,j)] / float(contactHistoNat1Count)
			if (i,j) in contactHistoNat2:
				norm_contactHistoNat2[(i,j)] = contactHistoNat2[(i,j)] / float(contactHistoNat2Count)
			if (i,j) in contactHistoNonNat:
				norm_contactHistoNonNat[(i,j)] = contactHistoNonNat[(i,j)] / float(contactHistoNonNatCount)
			if (i,j) in contactHistoFiltered:
				norm_contactHistoFiltered[(i,j)] = contactHistoFiltered[(i,j)] / float(contactHistoFilteredCount)
		
	
	assert writeContactHisto(length, norm_contactHistoAll, label + "_ALL_%s_contacts.dat"%outputdir)
	assert writeContactHisto(length, norm_contactHistoNat1, label + "_NATIVE1_%s_contacts.dat"%outputdir)
	assert writeContactHisto(length, norm_contactHistoNat2, label + "_NATIVE2_%s_contacts.dat"%outputdir)
	assert writeContactHisto(length, norm_contactHistoNonNat, label + "_NONNAT_%s_contacts.dat"%outputdir)
	assert writeContactHisto(length, norm_contactHistoFiltered, label + "_RES"+"-".join([str(i) for i in filterresidues])+"_%s_contacts.dat"%outputdir )
	
	print
	for key, value in sorted(norm_contactHistoFiltered.items(), key=lambda kv: kv[1], reverse=True):
		print key[0]+1, key[1]+1, value
	
	
	if plots:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as p
		import matplotlib.cm as cm
		from matplotlib.ticker import MultipleLocator, FormatStrFormatter
		
		plotContactMap(length, norm_contactHistoAll,      native1, native2, label, cutoff, chainDist, atoms, "ALL", label + "_ALL_%s_contacts.png"%outputdir)
		plotContactMap(length, norm_contactHistoNat1,     native1, native2, label, cutoff, chainDist, atoms, "NATIVE1", label + "_NATIVE1_%s_contacts.png"%outputdir)
		plotContactMap(length, norm_contactHistoNat2,     native1, native2, label, cutoff, chainDist, atoms, "NATIVE2", label + "_NATIVE2_%s_contacts.png"%outputdir)
		plotContactMap(length, norm_contactHistoNonNat,   native1, native2, label, cutoff, chainDist, atoms, "NONNAT", label + "_NONNAT_%s_contacts.png"%outputdir)
		plotContactMap(length, norm_contactHistoFiltered, native1, native2, label, cutoff, chainDist, atoms, "RES"+"-".join([str(i) for i in filterresidues]), label + "_RES"+"-".join([str(i) for i in filterresidues])+"_%s_contacts.png"%outputdir)
	
