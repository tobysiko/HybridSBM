import sys
from libprofsbm import splitModels

def run_splitPDB(arguments):
	if not type(arguments)==list: arguments = [arguments]
	inpdbname = arguments[0]
	pdblist = splitModels(inpdbname)
	return pdblist

if __name__=="__main__":
	print run_splitPDB(sys.argv[1:])