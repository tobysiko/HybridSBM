import sys, argparse

from libprofsbm import BB_atoms, splitModels, atomsFromPDB, getCutoffContacts, getShadowContacts, mixCutoffShadow, getResidueContacts, filterSuperEnsembleContacts, writeProfasiContactsFMULTIGAUSS, writeChimeraPseudobonds, plotMultiGaussianContacts, getLength, plotContactMap

#aromatic = ['PHE','TYR','TRP']


class parseSelectedRange(argparse.Action):

	def __call__(self, parser, namespace, values, option_string=None):
		print values
		return argparse.Action.__call__(self, parser, namespace, values, option_string)

	
	
	
	

def run_contactmap(sysargs):
	##################
	# HANDLE OPTIONS (-h for help)
	description="%prog prepares native contact files in XML format for PROFASI as modified by TS.\n Default values for options are given in parentheses."
	
	op = argparse.ArgumentParser(description=description)
	
	op.add_argument('inPDBFilename', help="structure file in PDB format. Comma-separate multiple files in order to use potential type FMULTIGAUSS and omit -S option.")
	op.add_argument('-c','--contactType', default="RESIDUE", choices=["CUTOFF", "SHADOW", "HYBRID", "RESIDUE"], help="type of contact definition: CUTOFF, SHADOW, HYBRID, RESIDUE  (%default)")
	op.add_argument('-a','--contactAtoms', default="CA", help="list of atom types for native contacts, e.g. CA,CB,C,N,O  (%default)")
	op.add_argument('-n','--chainSep', type=int, default=3, help="min number of residues between contacting atoms  (%default)")
	op.add_argument('-d','--distCutoff', type=float, default=6.0, help="distance cutoff for defining contacts (%default Ang)")
	op.add_argument('-m','--wellDepth', type=float, default=1.0, help="well depth at native distance (%default)")
	op.add_argument('-w','--wellWidth', type=float, default=0.5, help="GAUSS/DUALGAUSS: width of well (%default)")
	op.add_argument('-r','--repulDist', type=float, default=0.0, help="GAUSS/DUALGAUSS: excluded volume radius.  0=no repulsion; -x=radius is min distance per contact - x (FMULTIGAUSS)  (%default)")
	op.add_argument('-p','--repulSteep', type=float, default=1.0, help="GAUSS/DUALGAUSS: steepness of repulsion (%default)")
	op.add_argument('-e','--ensembleCutoff', type=float, default=1.0, help="when PDB file contains multiple models, retain contacts present in this fraction of models (%default)")
	op.add_argument('-N','--normEnergy', type=float, default=0.0, help="normalize contact energies to given amount per structure. Overrides options m and M. (%default=don't normalize)")
	op.add_argument('--gmxpath', default="/usr/local/gromacs/bin", help="path to Gromacs executables (SHADOW)")
	op.add_argument('--scmpath', default=".", help="path to SCM.jar  (SHADOW)")
	op.add_argument('-L','--label', default="", help="string that will be appended to the name 'SBM' in the Profasi energy output")
	op.add_argument('-S','--seqFrom', type=int, default=0, help="model number to select AA sequence for contacts. By default: take first model")
	op.add_argument('-x','--includeNonNative', action="store_true", help="include non-native contacts")
	op.add_argument('-R','--selectedRange', default="", help="specify residues ranges as strings,e.g. '1-20,29,30,42-56' add 'x' at the beginning to exclusively only consider these residues, otherwise contacts involving only one of these residues will be considered")
	op.add_argument('-M','--smoothGaussian',action='store_true',help='set function to epsilon if between min and max distance')
	op.add_argument('-q','--quiet',action='store_true',help='suppress all output to STDOUT')
	
	
	args = op.parse_args(sysargs)
	
	
	assert args.contactType in ["CUTOFF","SHADOW","HYBRID","RESIDUE"]
	assert args.distCutoff > 1.0, "cutoff in Angstrom!"
	
	restraintType = "sbm"
	
	contactPotential = "FMULTIGAUSS"
	if args.smoothGaussian:
		contactPotential = "FMULTIGSMOOTH"
	assert contactPotential in ["FMULTIGAUSS","FMULTIGSMOOTH"]
	
	selection = []
	selectExclusive = False
	if args.selectedRange!="":
		
		tmp = args.selectedRange.strip().split(",")
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
		if not args.quiet: print "Your selection of amino acid positions: ", args.selectedRange, selection
	
		
	consStatement = "_cons%s"%str(args.ensembleCutoff)
	
	args.contactAtoms = args.contactAtoms.split(",")
	
	if args.contactType == "SHADOW" and args.contactAtoms != ['CA']:
		print "SHADOW only allows CA atoms!"
		args.contactAtoms = ['CA']
	
	if args.label=="":
		args.label = "_"+str(args.inPDBFilename)[:-4]
	
	if not args.quiet: print "Atoms to consider for CUTOFF contacts:",args.contactAtoms
	
	assert all([ a in BB_atoms for a in args.contactAtoms])
	
	if args.contactType in ["HYBRID","SHADOW"]:
		if not args.quiet: print "Also considering CA atoms for SHADOW contacts."
	
	if args.normEnergy > 0 and not args.quiet: print "Normalizing contact energy per native structure to:", args.normEnergy
	
	assert args.inPDBFilename != "" and args.inPDBFilename != None
	
	atoms_models = []
	models_contact_lists = []
	models_scenergy_lists = []
	
	#atoms_models = extractModelAtoms(inPDBFilename)
	models = splitModels(args.inPDBFilename)
	
	if not args.quiet: print "\n\n%i model(s) found in %s."%(len(models), args.inPDBFilename)
	
	for m in models:
		con_c=[]
		con_s=[]
		
		atoms = atomsFromPDB(m)
		atoms_models.append(atoms)
		
		if not args.quiet: print "num atoms",len(atoms), m
		
		if args.contactType=="CUTOFF":
			m_contacts = getCutoffContacts(atoms, atomtypes=args.contactAtoms, cutoff=args.distCutoff, chainDist=args.chainSep, verbose=not args.quiet, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		elif args.contactType=="SHADOW":
			m_contacts = getShadowContacts(pdbfilename=m, outfilename=args.inPDBFilename[:-4]+".shadow", cutoff=args.distCutoff, chainDist=3, gmxpath=args.gmxpath, scmpath=args.scmpath, res2SS = None, selection=selection, selectExclusive=selectExclusive, verbose=not args.quiet)
		elif args.contactType=="HYBRID":
			con_c = getCutoffContacts(atoms, atomtypes=args.contactAtoms, cutoff=args.distCutoff, chainDist=args.chainSep, verbose=not args.quiet, res2SS = None, selection=selection, selectExclusive=selectExclusive)
			con_s = getShadowContacts(pdbfilename=m, outfilename=args.inPDBFilename[:-4]+".shadow", cutoff=args.distCutoff, chainDist=3, gmxpath=args.gmxpath, scmpath=args.scmpath, res2SS = None, selection=selection, selectExclusive=selectExclusive, verbose=not args.quiet)
			m_contacts = mixCutoffShadow(con_c, con_s)
		elif args.contactType=="RESIDUE":
			m_contacts = getResidueContacts(atoms, cutoff=args.distCutoff, chainDist=args.chainSep,verbose=not args.quiet, res2SS = None, selection=selection, selectExclusive=selectExclusive)
		else:
			print "unknown contact type: ",args.contactType
			sys.exit(1)
		
		models_contact_lists.append(m_contacts)
	
	if not args.quiet: print "%s ensemble contacts"%contactPotential
	finalcontacts = filterSuperEnsembleContacts(models_contact_lists, atoms_models, args.contactAtoms, args.ensembleCutoff, args.contactType,verbose=not args.quiet)
	
	for m in models:
		energies_m =[]
		models_scenergy_lists.append( energies_m  )
	
	pbond_legend = None
	
	outatoms = "".join([i for i in args.contactAtoms])
	
	width_label = "_w"+str(round(args.wellWidth,1))
	repul_label = "_R"+str(round(args.repulDist,1))
	nonnat_label = "_x"+str(args.includeNonNative)
	
	outfilename = "%s_%s_%s_c%s_n%s%s%s_%s%s%s_%s"%(args.inPDBFilename.split("/")[-1],contactPotential,args.contactType,str(args.distCutoff),str(args.chainSep),width_label,repul_label,outatoms,consStatement,nonnat_label,args.selectedRange)
	
	if args.contactType=="HYBRID" and not 'CA' in args.contactAtoms:
		args.contactAtoms.append('CA')
	
	if not args.quiet: print "Sequence from model",args.seqFrom
	
	outfilename = writeProfasiContactsFMULTIGAUSS(finalcontacts, atoms_models[args.seqFrom], args.contactAtoms, outfilename, args.repulDist, args.repulSteep, width=args.wellWidth, depth=args.wellDepth, norm=args.normEnergy,label=args.label,includeNonNative=args.includeNonNative,chainSep=args.chainSep,smoothGaussian=args.smoothGaussian, verbose=not args.quiet)
	
	writeChimeraPseudobonds(finalcontacts, outfilename[:-4]+".pseudobonds",legend=pbond_legend, verbose=not args.quiet)
	
	if not args.quiet: print "wrote Profasi XML contacts:",outfilename
	if not args.quiet: print "wrote Chimera pseudobonds:",outfilename[:-4]+".pseudobonds"
	
	if True:
		plotMultiGaussianContacts(finalcontacts,1.0,0.5,outfilename[:-4]+"_multiwell.png",smoothGaussian=args.smoothGaussian, verbose=not args.quiet)
		if not args.quiet: print "plotted all contact potentials to PNG:", outfilename[:-4]+"_multiwell.png"
	if True:
		nres = getLength(atoms)
		plotContactMap(nres, finalcontacts, args.inPDBFilename, args.distCutoff, args.chainSep, atoms, args.contactType, fname=outfilename[:-4]+".png", verbose=not args.quiet)
	
	return outfilename
	
if __name__=="__main__":
	print run_contactmap(sys.argv[1:])
