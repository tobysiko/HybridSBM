import argparse
from fetchPDB import run_fetchPDB
from cleanPDB import run_cleanPDB
from aggregatePDB import run_aggregatePDB
from contactmap import run_contactmap
from splitPDB import run_splitPDB

if __name__=="__main__":
	##################
	# HANDLE OPTIONS (-h for help)
	description = "%prog: runs automatic preparation from groups of PDB codes to SBM files ready for use in PROFASI."
	op = argparse.ArgumentParser(description=description)
	sbmstring=None
	op.add_argument('sbmstring', help="Provide PDB codes where s refers to the same sbm and p is a pdb code: e.g. s1p1,s1p2;s2p1,s2p2 produces two SBMs with two PDB files each")
	op.add_argument('-a','--allModels', action="store_true", dest="allmodels",help="When PDB file contains multiple models, use all of them for the consensus contact map. Default: only use first model.")
	op.add_argument('-l','--sbmlabels', help="Provide one label for each SBM in the form: label1,label2,label3 ...")
	op.add_argument('-e','--minetot', action="store_true", help="After regularization, use model with lowest Profasi energy. Default: use model with lowest RMSD from original.")
	op.add_argument('-q', '--quiet', action="store_true")
	args = op.parse_args()
	
	sbmstring_split1 = args.sbmstring.split(";")
	
	if args.sbmlabels: assert len(args.sbmlabels.split(",")) == len(sbmstring_split1)
	
	sbms = []
	for ai in xrange(len(sbmstring_split1)):
		a = sbmstring_split1[ai]
		codes = a.split(',')
		#print codes
		if args.sbmlabels==None:
			label = "consensus%i.pdb"%(ai)
		else:
			label = args.sbmlabels.split(",")[ai]
		
		if args.allmodels:
			pdblist = []
			for l in [run_splitPDB(run_cleanPDB(run_fetchPDB(c))) for c in codes]:
				pdblist.extend(l)
		else:
			pdblist = [run_splitPDB(run_cleanPDB(run_fetchPDB(c)))[0] for c in codes]
		print "pdblist:",pdblist
		
		agg_ar = [ label, ",".join(pdblist) ]
		
		if not args.quiet: agg_ar.append("-v")
		conspdb = run_aggregatePDB( agg_ar )
		
		cmap_args = [conspdb]
		if args.quiet: cmap_args.append("--quiet")
		cmap = run_contactmap(cmap_args)
		sbms.append(cmap)
		print cmap
