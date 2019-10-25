import sys, os, os.path, copy, math

from libprofsbm import *

if __name__=="__main__":
	print "takes a PDB file containing multiple models and each model with multiple chains and creates a new PDB file where each chain from each previous model is now a separate model with a single chain (default 'A')"
	
	inpdbname = sys.argv[1]
	outpdbname = sys.argv[2]
	
	print "splitting",inpdbname
	modelchains = splitModelsAndChains(inpdbname)
	print "joining",outpdbname
	joinModels(modelchains, outpdbname)
	
