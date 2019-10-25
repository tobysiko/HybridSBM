import sys,math,subprocess,os.path

def splitModels(pdbfilename):
	hasModels = False
	models = []
	model_nr = 0
	atomlines = []
	for line in open(pdbfilename).readlines():
		linetype = line[:6]
		
		if linetype == "MODEL ":
			hasModels = True
			model_nr = int(line.replace("MODEL",""))
			print "model nr ",model_nr
			atomlines = []
		elif linetype == "ENDMDL":
			assert model_nr != 0, model_nr
			modelfilename = pdbfilename[:-4] + "_m%i.pdb"%model_nr
			models.append(modelfilename)
			
			modelfile = open(modelfilename,"w")

			assert len(atomlines) != 0
			for l in atomlines:
				modelfile.write(l)
			modelfile.write("END   ")
			modelfile.close()
		elif linetype == "ATOM  ":
			atomlines.append(line)
		else:
			pass
	print len(models),"models found in",pdbfilename
	if hasModels:
		return models # returns a list of the new filenames, one model per file
	else:
		return [pdbfilename]

def joinModels(listOfFilenames, outfilename, cleanUp=True):
	print listOfFilenames
	outfile = open(outfilename, 'w')
	counter = 1
	for infilename in listOfFilenames:
		infile = open(infilename)

		spaces = (9-int(math.log10(counter))) * " "
		outfile.write("MODEL"+spaces+str(counter)+"\n")
		outfile.write("REMARK original_filename %s\n"%(infilename))
		for line in infile.readlines():

			if line[:4] == "ATOM":
				outfile.write(line)

		outfile.write("ENDMDL                                                                          \n")
		counter += 1
		infile.close()
	outfile.write("END                                                                             \n")
	outfile.close()
	if cleanUp:
		for l in listOfFilenames:
			 if os.path.exists(l): os.remove(l)

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


if __name__=="__main__":
	outname = sys.argv[1]
	doRegularization = sys.argv[2]
	if doRegularization in [1,"True","true","TRUE"]:
		doRegularization = True
	else:
		doRegularization = False
		
	if doRegularization: print "regularizing for Profasi!"
	listFiles = sys.argv[3:]
	
	splitFiles = []
	for f in listFiles:
		
		splitFiles.extend(splitModels(f))
	
	regFiles_rmsd = []
	regFiles_etot = []
	regFiles = []
	
	for f in splitFiles:
		print f
		if doRegularization:
			regularize(f,"~/workspace/PROFASI/app/bin/")
			print "adding:",f[:-4]+"_min_rmsd.pdb"
			regFiles_rmsd.append(f[:-4]+"_min_rmsd.pdb")
			regFiles_etot.append(f[:-4]+"_min_etot.pdb")
		else:
			regFiles.append(f)
		
	if doRegularization:
		joinModels(regFiles_rmsd, outname[:-4]+"_prf_min_rmsd.pdb")
		joinModels(regFiles_etot, outname[:-4]+"_prf_min_etot.pdb")
	else:
		joinModels(regFiles, outname)
	
	