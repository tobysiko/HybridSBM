#!/home/sikosek/Library/Enthought/Canopy_64bit/User/bin/pythonw
import sys,os,os.path,math,operator
import matplotlib.pyplot as p


def getBinIndex(bins, v):
	for i in xrange(len(bins)):
		if i==len(bins)-1:
			return i
		elif v >= bins[i] and v < bins[i+1]:
			#print v,bins[i],bins[i-1]
			return i

def d2r(degrees):
	if degrees != None:
		return degrees / (180/math.pi)
	else:
		return None

def r2d(radians):
	if radians != None:
		return radians * (180/math.pi)
	else:
		 return None
def tlc2olc(tlc):
	d = {"PHE":"F", "TYR":"Y", "TRP":"W", "HIS":"H", "ARG":"R", "LYS":"K","W6C":"W"}
	return d[tlc.upper()]
#assuming we are just outside the "pdb" folder that contains all the .pidat files


def correctF3D(dat, bins, binmidpoints, cdist_cut, F_scale, rep_F):
	for k in xrange(len(bins[2])):
		for j in xrange(len(bins[1])):	
		
			for i in xrange(len(bins[0])):
				if dat[(i,j,k)] != None:
					dat[(i,j,k)] /= F_scale
		
			for i in xrange(len(bins[0])):
				if binmidpoints[0][i] > cdist_cut:
					if dat[(i,j,k)] == None:
						dat[(i,j,k)] = 0.0
					elif dat[(i,j,k)] >= 0.0: 
						dat[(i,j,k)] = 0.0
		
			for i in xrange(len(bins[0])):
				if i>0 and i<len(bins[0])-1:
					if dat[(i,j,k)] == None and dat[(i-1,j,k)] != None and dat[(i+1,j,k)] != None:
						if dat[(i-1,j,k)] > rep_F: dat[(i-1,j,k)] = rep_F
						if dat[(i+1,j,k)] > rep_F: dat[(i+1,j,k)] = rep_F
						dat[(i,j,k)] = (dat[(i+1,j,k)] + dat[(i-1,j,k)]) / 2.0
		
			for i in xrange(len(bins[0])):
				if dat[(i,j,k)] == None and dat[(i+1,j,k)] == None and dat[(i+2,j,k)] == None:
					dat[(i,j,k)] = rep_F
			
		
			if dat[(0,j,k)]==None and dat[(1,j,k)] != None:
				if dat[(1,j,k)] > rep_F: dat[(1,j,k)]=rep_F
				diff = rep_F - dat[(1,j,k)] 
				dat[(0,j,k)] == dat[(1,j,k)] + diff/3.0
			elif dat[(0,j,k)]==None and dat[(1,j,k)] == None and dat[(2,j,k)] != None:
				if dat[(2,j,k)] > rep_F: dat[(2,j,k)]=rep_F
				diff = rep_F - dat[(2,j,k)]
				dat[(1,j,k)] == dat[(2,j,k)] + diff/3.0
				dat[(0,j,k)] == dat[(2,j,k)] + 2*diff/3.0
		
		
			for i in xrange(len(bins[0])-3):
				if dat[(i,j,k)] != None and dat[(i+1,j,k)] == None and dat[(i+2,j,k)] == None and dat[(i+3,j,k)] != None:
					if dat[(i,j,k)] > rep_F: dat[(i,j,k)] = rep_F
					if dat[(i+3,j,k)] > rep_F: dat[(i+3,j,k)] = rep_F
				
					diff = dat[(i,j,k)] - dat[(i+3,j,k)] 
					dat[(i+2,j,k)] = dat[(i+3,j,k)] + diff/3.0
					dat[(i+1,j,k)] = dat[(i+3,j,k)] + 2*diff/3.0
		
			for i in xrange(len(bins[0])):
				if dat[(i,j,k)] == None:
					dat[(i,j,k)] = rep_F
			
			for i in xrange(len(bins[0])):
				if i==0:
					if dat[(i,j,k)] != rep_F and dat[(i+1,j,k)] == rep_F:
						dat[(i,j,k)] = rep_F
				elif binmidpoints[0][i] < cdist_cut:
					if dat[(i,j,k)] != rep_F and dat[(i+1,j,k)] == rep_F and dat[(i-1,j,k)] == rep_F:
						dat[(i,j,k)] = rep_F

def correctF2D(dat, bins, binmidpoints, cdist_cut, F_scale, rep_F):
	for j in xrange(len(bins[1])):	
		
		for i in xrange(len(bins[0])):
			if dat[i][j] != None:
				dat[i][j] /= F_scale
		
		for i in xrange(len(bins[0])):
			if binmidpoints[0][i] > cdist_cut:
				if dat[i][j] == None:
					dat[i][j] = 0.0
				elif dat[i][j] >= 0.0: 
					dat[i][j] = 0.0
		
		for i in xrange(len(bins[0])):
			if i>0 and i<len(bins[0])-1:
				if dat[i][j] == None and dat[i-1][j] != None and dat[i+1][j] != None:
					if dat[i-1][j] > rep_F: dat[i-1][j] = rep_F
					if dat[i+1][j] > rep_F: dat[i+1][j] = rep_F
					dat[i][j] = (dat[i+1][j] + dat[i-1][j]) / 2.0
		
		for i in xrange(len(bins[0])):
			if dat[i][j] == None and dat[i+1][j] == None and dat[i+2][j] == None:
				dat[i][j] = rep_F
			
		
		if dat[0][j]==None and dat[1][j] != None:
			if dat[1][j] > rep_F: dat[1][j]=rep_F
			diff = rep_F - dat[1][j] 
			dat[0][j] == dat[1][j] + diff/3.0
		elif dat[0][j]==None and dat[1][j] == None and dat[2][j] != None:
			if dat[2][j] > rep_F: dat[2][j]=rep_F
			diff = rep_F - dat[2][j]
			dat[1][j] == dat[2][j] + diff/3.0
			dat[0][j] == dat[2][j] + 2*diff/3.0
		
		
		for i in xrange(len(bins[0])-3):
			if dat[i][j] != None and dat[i+1][j] == None and dat[i+2][j] == None and dat[i+3][j] != None:
				if dat[i][j] > rep_F: dat[i][j] = rep_F
				if dat[i+3][j] > rep_F: dat[i+3][j] = rep_F
				
				diff = dat[i][j] - dat[i+3][j] 
				dat[i+2][j] = dat[i+3][j] + diff/3.0
				dat[i+1][j] = dat[i+3][j] + 2*diff/3.0
		
		for i in xrange(len(bins[0])):
			if dat[i][j] == None:
				dat[i][j] = rep_F
		
		for i in xrange(len(bins[0])):
			if i==0:
				if dat[i][j] != rep_F and dat[i+1][j] == rep_F:
					dat[i][j] = rep_F
			elif binmidpoints[0][i] < cdist_cut:
				if dat[i][j] != rep_F and dat[i+1][j] == rep_F and dat[i-1][j] == rep_F:
					dat[i][j] = rep_F
		
		#for i in range(len(bins[0])-1,-1,-1):
		#	if binmidpoints[0][i] <= cdist_cut:
		#		if dat[i][j] != None:
		#			
		#			if dat[i][j] > rep_F:
		#				dat[i][j] = rep_F
		#			print i, dat[i][j]
		#		else:
		#			
		#			ii=i
		#			while dat[ii][j] == None and ii>0:
		#				ii -= 1
		#			
		#			
		#			ngaps = i-ii
		#			
		#			print i,ii,ngaps
		#			if i==0:
		#				
		#			if ngaps == 1:
		#				if dat[ii][j] > rep_F: dat[ii][j] = rep_F
		#				dat[i][j] = (dat[i+1][j] + dat[ii][j]) / 2.0
		#				print "ngap2",ii,dat[ii][j]
		#			elif ngaps==2:
		#				if dat[ii-1][j] > rep_F: dat[ii-1][j] = rep_F
		#				diff = dat[i][j] - dat[ii-1][j]
		#				dat[i-1][j] = dat[i][j]+ diff/3.0
		#				print i-1, dat[i-1][j]
		#				dat[ii][j] = dat[i-1][j]+ diff/3.0
		#				print ii,dat[ii][j]
		#			elif ngaps > 2:
		#				if dat[ii-1][j] == None: dat[ii-1][j] = rep_F
		#				diff = dat[i][j] - dat[ii-1][j]
		#				dat[i-1][j] = dat[i][j]+ diff/3.0
		#				print i-1, dat[i-1][j]
		#				dat[ii][j] = dat[i-1][j]+ diff/3.0
		#				print ii,dat[ii][j]
		#				for iii in range(ii-1,ii-1-ngaps-3,-1):
		#					dat[iii][j] = rep_F
		#					print iii,dat[iii][j]
							
					
					
				
					
				
				
				


def correctF2Dold(dat, bins, binmidpoints, cdist_cut, F_scale, rep_F):
	for j in xrange(len(bins[1])):	
		for i in xrange(len(bins[0])):
			#print i, bins[0][i]
			indices_scaled={}
			
			if binmidpoints[0][i] > cdist_cut:
				if j==0:print i,dat[i][j],
				#dat[i][j] = 0.8
				#continue
				if not i in indices_scaled: 
					indices_scaled[i]=True
					if dat[i][j] == None:
						dat[i][j] = 0.0
					elif dat[i][j] >= 0.0: 
						dat[i][j] = 0.0
					elif dat[i][j] < 0.0:
						
						dat[i][j] = dat[i][j] / F_scale
						
						
					else:
						print "???"
				if j==0: print dat[i][j]
				
				
				
			else:
				#print dat[i][j], dat[i+1][j], dat[i+2][j], dat[i+3][j], " --> ",
				#if dat[i][j] == None:
				if j==0:print i,dat[i][j],
				if dat[i][j] == None and dat[i+1][j] == None and dat[i+2][j] == None and dat[i+3][j] == None:
					if not i in indices_scaled: 
						indices_scaled[i]=True
						dat[i][j] = rep_F
			
				elif dat[i][j] == None and dat[i+1][j] == None and dat[i+2][j] == None and dat[i+3][j] != None:
					if not i+3 in indices_scaled: 
						indices_scaled[i+3]=True
						if dat[i+3][j] < 0.0:
							dat[i+3][j] = dat[i+3][j] / F_scale
						
						elif dat[i+3][j] > 0.0:
							dat[i+3][j] = max(dat[i+3][j] / F_scale, rep_F)
						else:
							print "i+3 = 0"
					
					if not i in indices_scaled: 
						indices_scaled[i]=True
						dat[i][j] = rep_F
					if not i+1 in indices_scaled: 
						indices_scaled[i+1]=True
						dat[i+1][j] = rep_F - abs(dat[i+3][j] - dat[i][j])/ 3.0
					if not i+2 in indices_scaled: 
						indices_scaled[i+2]=True
						dat[i+2][j] = rep_F - 2.0*abs(dat[i+3][j] - dat[i][j])/ 3.0
			
				elif dat[i][j] == None and dat[i+1][j] == None and dat[i+2][j] != None:
					if i>0:
						hival = dat[i-1][j]
					else:
						hival = rep_F
					
					if not i+2 in indices_scaled: 
						indices_scaled[i+2]=True
						if dat[i+2][j] < 0.0:
							dat[i+2][j] = dat[i+2][j] / F_scale
						elif dat[i+2][j] > 0.0:
							dat[i+2][j] = max(dat[i+2][j] / F_scale, rep_F)
						else:
							print "i+2 = 0"
					
					if not i in indices_scaled: 
						indices_scaled[i]=True
						dat[i][j] = hival - 2.0*abs(dat[i+2][j] - hival)/ 3.0
					
					if not i+1 in indices_scaled: 
						indices_scaled[i+1]=True
						dat[i+1][j]   = hival - abs(dat[i+2][j] - hival)/ 3.0
			
				elif dat[i][j] == None and dat[i+1][j] != None:
					
					if i>1:
						#assert i-2 in indices_scaled,i
						hival = dat[i-2][j]
					else:
						hival = rep_F
					
					if not i+1 in indices_scaled: 
						indices_scaled[i+1]=True
						if dat[i+1][j] < 0.0:
							dat[i+1][j] = dat[i+1][j] / F_scale
						elif dat[i+1][j] > 0.0:
							dat[i+1][j] = max(dat[i+1][j] / F_scale, rep_F)
						else:
							print "i+1 = 0"
					
					if not i in indices_scaled: 
						indices_scaled[i]=True
						if hival > dat[i+1][j]:
							dat[i][j] = hival - 2.0*abs(dat[i+1][j] - hival)/ 3.0
						else:
							dat[i][j] = hival + 2.0*abs(dat[i+1][j] - hival)/ 3.0
				
				elif dat[i][j] != None:
					if not i in indices_scaled: 
						indices_scaled[i]=True
						if dat[i][j] < 0.0:
							dat[i][j] = dat[i][j] / F_scale
						elif dat[i][j] > 0.0:
							if dat[i][j] / F_scale >= rep_F:
								dat[i][j] = rep_F
							else:
							
								dat[i][j] = dat[i][j] / F_scale
						else:
							print "i = 0"
				else:
					print "what else?"
				if j==0:print dat[i][j]
						
							
				#print dat[i][j], dat[i+1][j], dat[i+2][j], dat[i+3][j]



# EXAMPLE:
# python freeEnergyPiPi.py cmbi_pdbselect_90psi_rfact021_res20_20130126_n9796.txt PHE PHE cdist,theta,phi

pdblistfile = sys.argv[1]
res1 = sys.argv[2].upper()
res2 = sys.argv[3].upper()
columns = sys.argv[4].split(",")

if len(sys.argv)==6:
	givenpath=sys.argv[5]
else:
	givenpath=None

print columns
defaultBinSize = {"cdist":30,"theta":30,"phi":30,"psi":30,"zeta":30,"facc1":10,"facc2":10}

maxcdist = 12.0

cdist_cut = 5.0

rep_F = 1.0

F_scale = 30.3039723563

maxacc1 = None

fixedrange = True

showmat = False

ncontours = 30

pires = ["TRP","TYR","PHE","HIS","ARG","LYS","W6C"]

useW6C = True

if useW6C:
	w6clabel = "_W6C"
	if res1 == "TRP":
		res1 = "W6C"
		
	if res2=="TRP":
		res2 = "W6C"
else:
	w6clabel = ""

if res1 in pires[:4] and res2 in pires[:4]:
	maxTheta4Phi = 90.0
else:
	maxTheta4Phi = 90.0

#if len(sys.argv)==6:
#	nbins = [int(i) for i in sys.argv[5].split(",")]
#	assert len(nbins)== len(columns)
#else:
nbins = [defaultBinSize[c] for c in columns]

colnames = ["i","j","r1","r2","cdist","mdist","theta","phi","psi","zeta","facc1","facc2"]
assert all([c in colnames for c in columns])
assert "i" not in columns and "j" not in columns and "r1" not in columns and "r2" not in columns

colindices = []
for c in columns:
	colindices.append(colnames.index(c))

print colindices

assert len(colindices)==len(columns)


assert res1 in pires and res2 in pires

pidatlist = []

for line in open(pdblistfile).readlines():
	tmp = line.strip()
	
	if len(tmp.split())==1:
		if givenpath:
			fname = "%s/%s%s.pidat"%(givenpath,tmp[:-4],w6clabel)
		else:
			fname = "./pdb/%s%s.pidat"%(tmp,w6clabel)
		#print fname
		if os.path.exists(fname):
			pidatlist.append(fname)
		else:
			print "Error: could not find ",fname

rawdata = [[] for c in columns]
rawdata2 = [[] for c in columns]
#print len(rawdata)

print "reading .pidat files..."
count = 0
countNone = 0
countHiAngle=0
countLoAngle=0
ntotal = 0
datalabels=[]
for dat in pidatlist:
	sys.stdout.write("\r%s%%" % str(round(float(100.0*count)/len(pidatlist),2)))
	sys.stdout.flush()
	count+=1
	#print dat
	length = None
	labels = None
	
	tmp_orderpairs = []
	tmp_revpairs = []
	for line in open(dat).readlines():
		tmp = line.strip().split()
		#print tmp
		if len(tmp)==0:
			#print "?",line
			continue
		elif length==None:
			length = int(tmp[0])
		elif labels==None:
			labels = tmp
		else:
			assert len(tmp)==12,tmp
			
			checkOrderPair = (tmp[2]==res1 and tmp[3]==res2 and int(tmp[0]) < int(tmp[1]))
			checkReverseOrderPair = (tmp[2]==res2 and tmp[3]==res1 and int(tmp[0]) > int(tmp[1]))
			
			if checkOrderPair:
				tmp_orderpairs.append(tmp)
			elif checkReverseOrderPair:
				tmp_revpairs.append(tmp)
			#else:
			#	print "???",tmp
				
			#assert checkOrderPair != checkReverseOrderPair,tmp
			#if int(tmp[0]) < int(tmp[1]):
			#	tmp_orderpairs.append((tmp[0],tmp[1],tmp[2],tmp[3]))
			#elif int(tmp[0]) > int(tmp[1]):
			#	tmp_revpairs.append((tmp[0],tmp[1],tmp[2],tmp[3]))
			#	#if (tmp[1],tmp[0],tmp[3],tmp[2]) in tmp_orderpairs:
			#	#	tmp_orderpairs.remove((tmp[1],tmp[0],tmp[3],tmp[2]))
	
	#print tmp_orderpairs
	#print tmp_revpairs
	assert len(tmp_orderpairs) == len(tmp_revpairs),str(en(tmp_orderpairs))+" "+str(len(tmp_revpairs) )
	
	
	#sys.exit()
	
	for tmpi in tmp_orderpairs:
		
		for tmpj in tmp_revpairs:
			if tmpi[0]==tmpj[1] and tmpi[1]==tmpj[0]:
				#if  checkOrderPair or checkReverseOrderPair: 
				ntotal += 1
				#forward, tmpi
				if maxacc1 == None or (maxacc1 != None and tmpi[colnames.index("facc1")] != "None" and  float(tmpi[colnames.index("facc1")]) <= maxacc1): 
					#print float(tmp[colnames.index("facc1")])
					data_ok_i = True
					for c in colindices:
						if tmpi[c]=="None":
							data_ok_i = False
					#if data_ok_i and float(tmpi[4]) <= maxcdist:
					#	if checkOrderPair:datalabels.append(dat)
					#	if float(tmpi[colindices[1]]) <= 20.0: 
					#		#print dat,tmp
					#		countLoAngle += 1
					#	elif float(tmpi[colindices[1]]) >= 70.0: 
					#		#print dat,tmp
					#		countHiAngle += 1
					#	#for c in xrange(len(columns)):
					#	#	rawdata[c].append(float(tmpi[colindices[c]]))
					#		
					#elif not data_ok_j:
					#	
					#	countNone += 1
					#	continue
				#reverse, tmpj
				if maxacc1 == None or (maxacc1 != None and tmpj[colnames.index("facc1")] != "None" and  float(tmpj[colnames.index("facc1")]) <= maxacc1): 
					#print float(tmp[colnames.index("facc1")])
					data_ok_j = True
					for c in colindices:
						if tmpj[c]=="None":
							data_ok_j = False
					#if data_ok_j and float(tmpj[4]) <= maxcdist:
					#	if checkOrderPair:datalabels.append(dat)
					#	if float(tmpj[colindices[1]]) <= 20.0: 
					#		#print dat,tmp
					#		countLoAngle += 1
					#	elif float(tmpj[colindices[1]]) >= 70.0: 
					#		#print dat,tmp
					#		countHiAngle += 1
					#	#for c in xrange(len(columns)):
					#	#	rawdata2[c].append(float(tmpj[colindices[c]]))
					#elif not data_ok_j:
					#	countNone += 1
					#	continue
				
				# all good, proceed
				if data_ok_i and float(tmpi[4]) <= maxcdist and data_ok_j and float(tmpj[4]) <= maxcdist:
					datalabels.append(dat)
					if float(tmpi[colindices[1]]) <= 20.0: 
						#print dat,tmp
						countLoAngle += 1
					elif float(tmpi[colindices[1]]) >= 70.0: 
						#print dat,tmp
						countHiAngle += 1
					for c in xrange(len(columns)):
						rawdata[c].append(float(tmpi[colindices[c]]))
						rawdata2[c].append(float(tmpj[colindices[c]]))
				else:
					countNone += 1
				
	#for t in tmp_orderpairs
	#for c in xrange(len(columns)): 	assert len(rawdata[c])==len(rawdata2[c]), "%i %i"%(len(rawdata[c]),len(rawdata2[c]))
	#if len(tmp_orderpairs)>0: 
	#	print "unclaimed:",tmp_orderpairs
	#	for t in tmp_orderpairs:
	#		if t[2]==res1 and t[3]==res2 or t[2]==res2 and t[3]==res1:
	#			print "asymmetrical data !",dat
	#			
	#			#sys.exit()

print "lo:",countLoAngle
print "hi:",countHiAngle

ndata= len(rawdata[0])
ndata2 = len(rawdata2[0])

print len(rawdata),len(rawdata2), ndata, ndata2,len(datalabels),ntotal,ntotal-ndata

assert ndata==len(datalabels)==ndata2



mins = [10e10 for c in columns]
maxs = [-10e10 for c in columns]

if fixedrange:
	for c in xrange(len(columns)):
		if columns[c] == "cdist":
			mins[c] = 3.0
			maxs[c] = maxcdist
		elif columns[c] in ["theta","phi"]:
			mins[c] = 0.0
			maxs[c] = 90.0
		elif columns[c] in ["psi","zeta"]:
			mins[c] = 0.0
			maxs[c] = 180.0
else:

	
	#print rawdata
	for r in xrange(ndata):
	
		for c in xrange(len(columns)):
			if rawdata[c][r] < mins[c]:
				mins[c] = rawdata[c][r]
				#print rawdata[c][r]
			if rawdata[c][r] > maxs[c]:
				maxs[c] = rawdata[c][r]

print "mins",mins
print "maxs",maxs



binsizes = [(maxs[c]-mins[c])/nbins[c] for c in xrange(len(columns))]
bins = [[round(mins[c] + n*binsizes[c],6) for n in xrange(nbins[c])] for c in xrange(len(columns))]
binmidpoints = [[round(bins[c][n] + (binsizes[c]/2.0),6) for n in xrange(nbins[c])] for c in xrange(len(columns))]


print "nbins",nbins
print "binsizes",binsizes
print "bins",bins,len(bins[0]),len(bins[1])
print "binmidpoints",binmidpoints

histogram = {}
histogram2 = {}
nskip = 0
for r in xrange(ndata):
	key = []
	key2=[]
	skip=False
	
	for c in xrange(len(columns)):
		val = rawdata[c][r]
		if val < mins[c]:# or val > maxs[c]:
			print "skip:",val
			skip = True
		b = bins[c] 
		bind = getBinIndex(b,val)
		key.append(bind)
		#if c==0 and bind == 0: print c,val,bind
		
		val2 = rawdata2[c][r]
		if val2 < mins[c]:# or val > maxs[c]:
			print "skip:",val2
			skip = True
		b = bins[c] 
		bind = getBinIndex(b,val2)
		key2.append(bind)
		#if val!=val2: print val,val2
	if skip: 
		nskip += 1
		continue
	key = tuple(key)
	key2=tuple(key2)
	#print key,bins[0][key[0]],bins[1][key[1]]
	if key in histogram:
		histogram[key] += 1
	else:
		histogram[key] = 1
	
	if key2 in histogram2:
		histogram2[key2] += 1
	else:
		histogram2[key2] = 1
		

print "data outside fixed range:",nskip

ndata -= nskip
#for h in sorted(histogram.keys()):
#	print h,histogram[h],bins[0][h[0]],bins[1][h[1]]

#RT = 8.3144621 * 300 #J / mol
RT = 1.9872041 * 10E-3 * 300 # kcal / mol
dat = []
Pmax = -1.0


if len(columns)==2:


	for i in xrange(len(bins[0])):
		tmp = []
		for j in xrange(len(bins[1])):
		
			if (i,j) in histogram:
				d=float(ndata)#1.0#bins[0][i]*bins[0][i]*math.sin( d2r(bins[1][j]) )
				P=histogram[(i,j)]/d
				if P > Pmax: Pmax = P
				tmp.append(P)
			else:
				tmp.append(0.0)
		dat.append(tmp)
	
	dat = [[dat[i][j]/Pmax if dat[i][j]!=0.0 else 0.0 for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
	dat = [[-RT*math.log(dat[i][j]/float(bins[0][i]*bins[0][i]*math.sin( d2r(bins[1][j]) ) )) if dat[i][j]!=0.0 else 0.0 for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
	maxdat = max([max(i)for i in dat])

	dat = [[ 1.0 - (dat[i][j]/maxdat) if dat[i][j]!=0.0 else 0.0 for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
	assert len(rawdata[0])==len(rawdata[1])

	p.contourf(binmidpoints[1], binmidpoints[0],dat , ncontours)
	p.title("%s-%s (n=%i; x=%i)"%(res1,res2,ndata,countNone))
	p.xlabel(columns[1])
	p.ylabel(columns[0])
	p.xlim(mins[1],maxs[1])
	p.ylim(mins[0],maxs[0])
	cbar = p.colorbar()
	cbar.set_label("1.0-(Fij/Fmax); Fij=-RTln(Pij/Pmax); Pij=nij/(r^2sin(theta))")
	p.show()
elif len(columns)==3:
	
	probs = {}
	for i in xrange(len(bins[0])): 
		
		for j in xrange(len(bins[1])):
			
			if binmidpoints[1][j] <= maxTheta4Phi:
			
				for k in xrange(len(bins[2])):
					if (i,j,k) in histogram:
						n = histogram[(i,j,k)]
					else:
						n = 0
				
					if (i,j,k) in histogram2:
						n2 = histogram2[(i,j,k)]
					else:
						n2 = 0
				
					if n==0 and n2==0:
						P=0.0
					else:
						g = n / (binmidpoints[0][i]*binmidpoints[0][i] * math.sin(d2r(binmidpoints[1][j])) * math.sin(d2r(binmidpoints[2][k])) )
						g2 = n2 / (binmidpoints[0][i]*binmidpoints[0][i] * math.sin(d2r(binmidpoints[1][j])) * math.sin(d2r(binmidpoints[2][k])) )
						P = g+g2 #* math.sin(d2r(bins[2][k])) * binsizes[2]
				
					probs[(i,j,k)] = P
	
	Pmax = max(probs.iteritems(), key=operator.itemgetter(1))[1]
	#print Pmax
	for l in probs.keys():
		#print probs[l]
		if probs[l] != 0.0:
			#probs[l] = -RT * math.log( probs[l] )
			probs[l] = -RT * math.log( probs[l]/Pmax )
		else:
			probs[l] = None
	
	#correctF3D(probs, bins, binmidpoints, cdist_cut, F_scale, rep_F)
	
	Fmax = max(probs.iteritems(), key=operator.itemgetter(1))[1]
	
	prf_file = open("%s%s_r%s_%s_%s_t%s_%s_%s_p%s_%s_%s_%s%s.pidat"%(pdblistfile.split(".")[0],w6clabel,str(mins[0]),str(maxs[0]),str(nbins[0]),str(mins[1]),str(maxs[1]),str(nbins[1]),str(mins[2]),str(maxs[2]),str(nbins[2]), tlc2olc(res1), tlc2olc(res2)), "w")
	
	for l in sorted(probs.keys()):
		i=l[0]
		j=l[1]
		k=l[2]
		
		if probs[l] != None:
			probs[l] = 1.0 - probs[l]/Fmax
		else:
			probs[l] = 0.0
		
		prf_file.write("\t".join([ str(bins[0][i]), str(bins[1][j]), str(bins[2][k]), str(probs[l])  ]) + "\n")
	
	print "Pmax,Fmax:",Pmax,Fmax
	
	prf_file.close()
	
	r_a = []
	th_a = []
	ph_a = []
	s_a=[]
	
	r_r = []
	th_r = []
	ph_r = []
	s_r=[]
	
	for l in sorted(probs.keys()):
		i=l[0]
		j=l[1]
		k=l[2]
		v=probs[l]
		#print l
		#if v < 0.5:
		r_r.append(i)
		th_r.append(j)
		ph_r.append(k)
		s_r.append(v)
		#elif v >= 0.5:
		#	r_a.append(i)
		#	th_a.append(j)
		#	ph_a.append(k)
		#	s_a.append(v)
	
	if False:
		from enthought.mayavi import mlab
		fig1 = mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0),size=(800,800))
		#mlab.points3d(ph_a, th_a, r_a, s_a,opacity=0.1,scale_factor=3)#, color=(0,0,1))
		mlab.points3d(ph_r, th_r, r_r, s_r,opacity=0.3,scale_factor=2)#, color=(1,0,0))
		mlab.axes(ranges=[0.0, 90.0, 0.0, 90.0, 3.0, 9.0])
		#mlab.orientation_axes()
		mlab.xlabel("phi")
		mlab.ylabel("theta")
		mlab.zlabel("centroid\ndistance")
		#mlab.colorbar(nb_labels=2,nb_colors=2,label_fmt='')
		mlab.colorbar()
		#mlab.text(0.1,0.1,"attraction",color=(0,0,1),width=0.1)
		#mlab.text(0.8,0.1,"repulsion",color=(1,0,0),width=0.1)
		mlab.text(0.1,0.8,"%s-%s (n=%i; x=%i); (%s)"%(res1,res2,ndata,countNone,pdblistfile),width=0.8)
		mlab.title("Free Energy (%s,%s)"%(res1,res2),size=0.3,height=0.7,figure=fig1)
		viewdist = 120
		elevation = 60  # angle or 'rotate'
		azimuth = 180+45 # angle or 'rotate'
		mlab.view(distance=viewdist,elevation=elevation, azimuth=azimuth)
	
		mlab.savefig(pdblistfile+"_%s_%s_r_theta_phi.png"%(res1,res2))
	
		mlab.show()
	
		#sys.exit()
	##############
	
	for i in xrange(len(bins[0])):
		tmp = []
		for j in xrange(len(bins[1])):
			sum_P = 0.0
			for k in xrange(len(bins[2])):
				if (i,j,k) in histogram:
					n = histogram[(i,j,k)]
				else:
					n = 0
				
				if (i,j,k) in histogram2:
					n2 = histogram2[(i,j,k)]
				else:
					n2 = 0
				
				if n==0 and n2==0:
					P=0.0
				else:
					g = n / (binmidpoints[0][i]*binmidpoints[0][i] * math.sin(d2r(binmidpoints[1][j])) * math.sin(d2r(binmidpoints[2][k])) )
					g2 = n2 / (binmidpoints[0][i]*binmidpoints[0][i] * math.sin(d2r(binmidpoints[1][j])) * math.sin(d2r(binmidpoints[2][k])) )
					P = (g+g2) * math.sin(d2r(bins[2][k])) * binsizes[2]
				
				sum_P += P
			
			tmp.append(sum_P)
			
		dat.append(tmp)
	
	Pmax = max([max(i)for i in dat])
	
	dat = [[dat[i][j]/Pmax for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
	
	dat = [[ -RT*math.log(dat[i][j]) if dat[i][j]!=0.0 else None for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
	
	Fmax = max([max(i)for i in dat])
	Fmin = min([min([j for j in i if j != None])for i in dat if not [j for j in i if j != None] == []])
	print "Theta fmax fmin:",Fmax,Fmin
	
	
	#correctF2D(dat, bins, binmidpoints, cdist_cut, F_scale, rep_F)
	
	#Fmax = max([max(i)for i in dat])
	#Fmin = min([min(i)for i in dat])
	
	dat = [[ 1.0 - dat[i][j]/Fmax if dat[i][j]!=None else 0.0 for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
	#print "Pmax,Fmax:",Pmax,Fmax
	
	outfile = open(pdblistfile+"_%s_%s_r_theta.dat"%(res1,res2),"w")
	outfile.write("d=r,l=th\t" + "\t".join([str(dd) for dd in binmidpoints[1]])+"\n")
	for d in xrange(len(dat)):
		outfile.write(str(binmidpoints[0][d])+"\t"+"\t".join([str(dd) for dd in dat[d]])+"\n")
	outfile.close()
	
	#####
	dat2 = []
	
	for i in xrange(len(bins[0])):
		tmp = []
		for j in xrange(len(bins[2])):
			sum_P = 0.0
			for k in xrange(len(bins[1])):
				if (i,k,j) in histogram:
					n=histogram[(i,k,j)]
				else:
					n=0
				if (i,k,j) in histogram2:
					n2=histogram2[(i,k,j)]
				else:
					n2=0
					
				if binmidpoints[1][k] <= maxTheta4Phi: # theta restricted to below 15 degrees!
					
					if n==0 and n2==0:
						P=0.0
					else:	
						g = n / (binmidpoints[0][i]*binmidpoints[0][i] * math.sin(d2r(binmidpoints[2][j])) * math.sin(d2r(binmidpoints[1][k])) )
						g2 = n2 / (binmidpoints[0][i]*binmidpoints[0][i] * math.sin(d2r(binmidpoints[2][j])) * math.sin(d2r(binmidpoints[1][k])) )
						#if g!=g2:print g,g2
						P = (g+g2) * math.sin(d2r(bins[1][k])) * binsizes[1]
				
				sum_P += P
			
			tmp.append(sum_P)
			
		dat2.append(tmp)
	
	Pmax2 = max([max(i)for i in dat2])
	
	dat2 = [[dat2[i][j]/Pmax2  for j in xrange(len(bins[2]))] for i in xrange(len(bins[0]))]
	
	dat2 = [[ -RT*math.log(dat2[i][j]) if dat2[i][j]!=0.0 else None for j in xrange(len(bins[2]))] for i in xrange(len(bins[0]))]
	
	Fmax2 = max([max(i)for i in dat2])
	Fmin2 = min([min([j for j in i if j != None])for i in dat2 if not [j for j in i if j != None] == []])
	print "Phi fmax fmin:",Fmax2,Fmin2
	
	#correctF2D(dat2, bins, binmidpoints, cdist_cut, F_scale, rep_F)
	
	#Fmax2 = max([max(i)for i in dat])
	#Fmin2 = min([min(i)for i in dat])
	
	
	dat2 = [[ 1.0 - dat2[i][j]/Fmax2 if dat2[i][j]!=None else 0.0 for j in xrange(len(bins[2]))] for i in xrange(len(bins[0]))]
	print "Pmax,Fmax:",Pmax2,Fmax2
	
	#dat2 = [[ 1.0 - dat2[i][j]/Fmax if dat2[i][j]!=None else 0.0 for j in xrange(len(bins[2]))] for i in xrange(len(bins[0]))]
	
	
	outfile = open(pdblistfile+"_%s_%s_r_phi.dat"%(res1,res2),"w")
	outfile.write("d=r,l=phi\t" + "\t".join([str(dd) for dd in binmidpoints[2]])+"\n")
	for d in xrange(len(dat2)):
		outfile.write(str(binmidpoints[0][d])+"\t"+"\t".join([str(dd) for dd in dat2[d]])+"\n")
	outfile.close()
	
	#####
	
	#fig = p.Figure(figsize=(14, 6),dpi=200)
	fig, axes = p.subplots(ncols=2,figsize=(16, 6),dpi=200)
	#p.setp(axes.flat, aspect="auto", adjustable='box-forced')
	#p.subplots_adjust(wspace=1)
	fig.suptitle("%s-%s (n=%i; x=%i); (%s)"%(res1,res2,ndata,countNone,pdblistfile))
	#ax1 = p.subplot(1, 2, 1,adjustable='box-forced',aspect="equal",axes.flat)
	axes[0].set_title("w(r,theta);rep=%s;cdist_cut=%s;Fscale=%s"%(str(rep_F),str(cdist_cut),str(F_scale)) )
	#p.autoscale(enable=True, axis=u'both', tight=True)
	#p.tight_layout(pad=2.08)
	#ax1.set_aspect("equal")
	if showmat:
		#dat = [[ dat[i][j] if dat[i][j]!=None else Fmax for j in xrange(len(bins[1]))] for i in xrange(len(bins[0]))]
		cont = axes[0].matshow(dat,origin='lower',aspect='equal')
	else:
		cont = axes[0].contourf(binmidpoints[1], binmidpoints[0],dat , ncontours, origin='lower',extent=(0,1,0,1))
	axes[0].set_xlabel(columns[1])
	axes[0].set_ylabel(columns[0])
	if not showmat:
		axes[0].set_xlim(mins[1],maxs[1])
		axes[0].set_ylim(mins[0],maxs[0])
	
	cbar = p.colorbar(cont,ax=axes[0])
	#cbar.set_label("")
	#cont.set_clim(Fmin,Fmax)
	
	
	
	
	#ax2 = p.subplot(1, 2, 2)
	axes[1].set_title("w(r,phi)theta<=%s;rep=%s;cdist_cut=%s;Fscale=%s"%(str(maxTheta4Phi),str(rep_F),str(cdist_cut),str(F_scale)))
	if showmat:
		#dat2 = [[ dat2[i][j] if dat2[i][j]!=None else Fmax for j in xrange(len(bins[2]))] for i in xrange(len(bins[0]))]
		cont2 = axes[1].matshow(dat2,origin='lower',aspect='equal')
	else:
		cont2 = axes[1].contourf(binmidpoints[2], binmidpoints[0],dat2 , ncontours)
	#p.title("%s-%s (n=%i; x=%i)"%(res1,res2,ndata,countNone))
	axes[1].set_xlabel(columns[2])
	axes[1].set_ylabel(columns[0])
	
	if not showmat:
		axes[1].set_xlim(mins[2],maxs[2])
		axes[1].set_ylim(mins[0],maxs[0])
	cbar2 = p.colorbar(cont2,ax=axes[1])
	cbar2.set_label("")
	#cont2.set_clim(Fmin2,Fmax2)
	
	savename = "%s%s_r%s_%s_%s_t%s_%s_%s_p%s_%s_%s_%s%s.png"%(pdblistfile.split(".")[0],w6clabel,str(mins[0]),str(maxs[0]),str(nbins[0]),str(mins[1]),str(maxs[1]),str(nbins[1]),str(mins[2]),str(maxs[2]),str(nbins[2]), tlc2olc(res1), tlc2olc(res2))
	p.savefig(savename)#pdblistfile+".png"%(res1,res2))
	#p.show()
	
	
