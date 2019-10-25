import sys, os, re

from libprofsbm import joinModels, getTemps


def run_exTraj(loc, outfolder):
	basename = os.path.split(loc)[1]
	
	# number of replicas (trajectories) in this simulation folder
	counter = 0
	for d in os.listdir(loc):
		m = re.match(r'n\d', d)
		
		if m != None: 
			counter += 1
	n = counter 
	
	print n,"simulation folders found"
	
	with open( os.path.join(loc,"n0/rtkey")) as f:
		print f.readline()#discard
		tmphead = [c.strip().replace(" ","_") for c in f.readlines()]
		colhead = ",".join( [tmphead[0]] + ["replica"] + tmphead[1:] )
		print colhead
	
	temps = getTemps(os.path.abspath(os.path.join(loc, "n0", "temperature.info")))
	minenpdb_list = []
	currentpdb_list = []
	outfilename = os.path.join(outfolder, "%s_n%i.csv"%(basename, n))
	with open(outfilename, "w") as outf:
		outf.write(colhead+"\n")
		for i in xrange(n):
			subf = "n%i"%i
			full = os.path.abspath(os.path.join(loc, subf))
			assert os.path.exists(full)
			with open(os.path.join(full,"rt")) as inpf:
				lines = inpf.readlines()
				
				for line in lines:
					linesplit = line.strip().split()
					
					if len(linesplit) > 2:
						assert len(linesplit) == len(tmphead), "%i,%i"%(len(linesplit), len(tmphead))
						newline = ""
						newline = ",".join( [linesplit[0]] + [str(i)] + [ str(temps[int(linesplit[1])][0]) ] + linesplit[2:])
						if line==lines[-1]: print newline
						outf.write(newline+"\n")
			minen = "%s/n%i/minen.pdb"%(os.path.abspath(loc),i)
			current = "%s/n%i/current.pdb"%(os.path.abspath(loc),i)
			if os.path.exists(minen):
				minenpdb_list.append(minen)
			else:
				print minen, "not found!"
			if os.path.exists(current):
				currentpdb_list.append(current)
			else:
				print current, "not found!"
	
	outname = "%s/%s_n%i_minen.pdb"%(os.path.join(outfolder), basename,n)
	joinModels(minenpdb_list, outname)
	print "Wrote minen PDB models to",outname
	
	outname = "%s/%s_n%i_current.pdb"%(os.path.join(outfolder), basename,n)
	joinModels(currentpdb_list, outname)
	print "Wrote current PDB models to",outname

if __name__=="__main__":
	if len(sys.argv)>=2:
		loc = os.path.abspath(sys.argv[1])
		print loc
	else:
		loc = os.getcwd()
	
	if len(sys.argv)==3:
		outfolder = sys.argv[2]
	else:
		outfolder = "."
	
	run_exTraj(loc, outfolder)
	