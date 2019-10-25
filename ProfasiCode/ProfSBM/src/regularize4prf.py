import sys
from libprofsbm import regularize

def run_regularize4prf(arguments):
	inpdb = arguments[1]
	prf_loc="~/workspace/PROFASI/app/bin/"
	print "Looking for "+prf_loc+"prf_convert ..."
	regularize(inpdb, prf_loc)

if __name__=="__main__":
	run_regularize4prf(sys.argv)
