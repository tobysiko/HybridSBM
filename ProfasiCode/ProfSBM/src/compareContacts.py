import sys,ast,math,numpy
import operator

def plotContactMap(n, contacts, title):
	import matplotlib.pyplot as p
	import matplotlib.cm as cm
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	
	minP = min([contacts[i] for i in contacts.keys()])
	maxP = max([contacts[i] for i in contacts.keys()])
	
	matrix = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	#matrix1 = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	#matrix2 = [[0.0 for j in xrange(n)] for i in xrange(n) ]
	
	for i in range(1,n+1):
		for j in range(1,n+1):
			if (i,j) in contacts:
				matrix[i][j] = contacts[(i,j)]
			
	
	#title = "%s; cutoff=%s; csep=%i; n=%i\n%s"%(inPDBFilename,str(round(cutoff,2)),chainSep,len(contacts),contactType)
	fig = p.figure(figsize=(5,5))
	ax = p.subplot(111)
	imax = p.matshow(matrix, cmap=cm.jet, aspect='equal')
	majorLocator   = MultipleLocator(10)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator   = MultipleLocator(1)
	ax.xaxis.set_major_locator(majorLocator)
	ax.xaxis.set_major_formatter(majorFormatter)
	ax.xaxis.set_minor_locator(minorLocator)
	ax.yaxis.set_major_locator(majorLocator)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.yaxis.set_minor_locator(minorLocator)
	p.title(title)
	cbar = p.colorbar()
	p.xlabel("residue position")
	p.ylabel("residue position")
	cbar.set_label("")
	#imax.get_axes().set_axisbelow(True)
	#p.savefig(figsavename)
	p.show()
def writeContactHisto(n, histo, filename):
	f = open(filename,"w")
	histostring = "{"
	#histo.values().sort(key=lambda x: x[col], reverse=True)
	for i,v in sorted(histo.items(), key=operator.itemgetter(1),reverse=True):
		histostring += "%s:%s,\n"%(i,histo[i])
	histostring += "}"
	print histostring
	f.write(str(n) + "\n")
	f.write(str(histostring) + "\n")
	f.close()
	return True

# module load intel/15.0 python/2.7.8




if len(sys.argv) > 2:
	f1 = sys.argv[1]
	f2 = sys.argv[2]
	outname = sys.argv[3]

	file1 = open(f1)
	file2 = open(f2)

	n = int(file1.readline())
	n2 = int(file2.readline())
	assert n==n2

	c1 = ast.literal_eval(file1.readline().strip())
	c2 = ast.literal_eval(file2.readline().strip())

	correlation = {}

	PP = []
	Pi = []
	Pj = []
	for i in range(1,n+1):
		for j in range(1,n+1):
			if i<j:
				if (i,j) in c1:
					P1 = c1[(i,j)]
				else:
					P1 = 0.0
			
				if (i,j) in c2:
					P2 = c2[(i,j)]
				else:
					P2 = 0.0
			
				if P1 != 0.0 and P2 != 0.0:
					#PP.append(P1*P2)
					Pi.append(P1)
					Pj.append(P2)
					correlation[(i,j)] = P1*P2 / math.pow(P1,2)
	
	avgP1=numpy.mean(Pi)
	avgP2=numpy.mean(Pj)
	
	if False:
		maxCorr = max([correlation[i] for i in correlation])
	for i in xrange(n):
		for j in xrange(n):
			if (i,j) in correlation:
				PP.append((P1-avgP1)*(P2-avgP2))
				#correlation[(i,j)] = (2*correlation[(i,j)]/maxCorr)-1
	
	print "matrix correlation E((Pi-<Pi>)(Pj-<Pj>)) / (s(Pi)*s(Pj)) =", numpy.mean(PP) / (numpy.std(Pi)*numpy.std(Pj))
	
	print "contacts PiPj/Pi^2"
	
	assert writeContactHisto(n, correlation, outname)
else:
	f1 = sys.argv[1]
	file1 = open(f1)
	n = int(file1.readline())
	correlation = ast.literal_eval(file1.readline().strip())


plotContactMap(n, correlation, f1)
		
