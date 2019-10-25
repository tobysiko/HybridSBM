import sys
from libprofsbm import downloadPDB

def run_fetchPDB(arguments):
	if not type(arguments)==list: arguments = [arguments]
	pdbcode = arguments[0]
	#print pdbcode
	if not downloadPDB(pdbcode): 
		print "download failed!"
		return None
	else:
		return pdbcode.upper()+".pdb"

if __name__=="__main__":
	print run_fetchPDB(sys.argv[1:])
