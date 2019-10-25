import sys, argparse
from libprofsbm import splitModels, joinModels, regularize


def run_aggregatePDB(arguments):
	description="%prog aggregates multiple PDB files into one."
	op = argparse.ArgumentParser(description=description)
	
	op.add_argument('outname', help="")
	op.add_argument('doRegularization', action="store_true", help="")
	op.add_argument('listFiles', help="")
	op.add_argument('--minetot', action="store_true", help="")
	op.add_argument('-v','--verbose',help="")
	
	args = op.parse_args(arguments)
	
	
		
	if args.doRegularization and args.verbose: print "regularizing for Profasi!"
	
	splitFiles = []
	
	for f in args.listFiles.split(","):
		
		splitFiles.extend(splitModels(f))
	
	regFiles_rmsd = []
	regFiles_etot = []
	regFiles = []
	
	for f in splitFiles:
		if args.verbose: print f
		if args.doRegularization:
			regularize(f,"~/workspace/PROFASI/app/bin/", options="-ll 0", verbose=args.verbose)
			
			
			if args.minetot:
				if args.verbose: print "adding:",f[:-4]+"_min_rmsd.pdb"
				regFiles_etot.append(f[:-4]+"_min_etot.pdb")
			else:
				if args.verbose: print "adding:",f[:-4]+"_min_rmsd.pdb"
				regFiles_rmsd.append(f[:-4]+"_min_rmsd.pdb")
		else:
			regFiles.append(f)
	
	if args.doRegularization:
		if args.minetot:
			joinModels(regFiles_etot, args.outname[:-4]+"_prf_min_etot.pdb")
			return args.outname[:-4]+"_prf_min_etot.pdb"
		else:
			joinModels(regFiles_rmsd, args.outname[:-4]+"_prf_min_rmsd.pdb")
			return args.outname[:-4]+"_prf_min_rmsd.pdb"
	else:
		joinModels(regFiles, args.outname)
		return args.outname
	

if __name__=="__main__":
	print run_aggregatePDB(sys.argv[1:])