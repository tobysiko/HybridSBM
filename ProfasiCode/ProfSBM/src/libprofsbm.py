from decimal import Decimal
#from matplotlib.collections import PatchCollection
from random import uniform
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.cm as cm
import matplotlib.pyplot as p

import numpy.linalg

import sys, glob, os, math, shutil, copy, subprocess, random
import numpy as np
import pandas as pd

# UNITS

# Temperature conversion factor: 1 profasi units in Kelvin.
pru_in_kelvin = 2000.0 / 3
# Temperature conversion factor: 1 Kelvin in profasi units
kelvin_in_pru = 1.0 / pru_in_kelvin
# Energy conversion factor: 1 profasi unit in kcal/mol
prf_energy_in_kcal_per_mol = 1.9858 * pru_in_kelvin / 1000
# Energy conversion factor: 1 kcal/mol in profasi unit
kcal_per_mol_in_prf_energy = 1.0 / prf_energy_in_kcal_per_mol  # added by TS
k_B = 0.0019872041
kBprf = k_B / kcal_per_mol_in_prf_energy * kelvin_in_pru

BB_atoms = ["C", "CA", "O", "N", "OXT", "CB"]

AA_3_to_1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
AA_1_to_3 = {'G':'GLY','A':'ALA','V':'VAL','L':'LEU','I':'ILE','C':'CYS','M':'MET','F':'PHE','Y':'TYR','W':'TRP','P':'PRO','S':'SER','T':'THR','N':'ASN','Q':'GLN','D':'ASP','E':'GLU','H':'HIS','K':'LYS','R':'ARG'}

def getTemps(filename):
	t = open(filename)
	h = t.readline().strip("#").strip().split()
	d = {}
	for line in t.readlines():
		tmp = line.strip().split()
		index = int(tmp[0])
		Tprf = float(tmp[1])
		Tkel = float(tmp[2])
		beta = float(tmp[3])
		d[index] = (Tprf, Tkel, beta)
	t.close()
	return d
def etotPerTemp(df, temps):
	etotPerTemp = {} 
	for t in temps:
		etotPerTemp[t] = []
	
	for i in xrange(len(df["temperature"])):
		t = df["temperature"][i]
		e = df["Etot"][i]
		etotPerTemp[t].append(e)
	return etotPerTemp

def getBinIndex(bins, v):
	for i in xrange(len(bins)):
		if i == len(bins) - 1:
			return i
		elif v >= bins[i] and v < bins[i + 1]:
			return i

def punb(histoPerTemp3D, f, ct0, countPerTemp, nEtotbins, nQAbins, nQBbins):
	p_unb = [[ [0.0 for b3 in xrange(nEtotbins)] for b2 in xrange(nQBbins)] for b1 in xrange(nQAbins)]
	
	for q1 in xrange(nQAbins):
		for q2 in xrange(nQBbins):
			for e in xrange(nEtotbins):
				enu = sum([histoPerTemp3D[t][q1][q2][e] for t in xrange(len(countPerTemp))])
				den = sum([countPerTemp[t] * f[t] * ct0[e][t] for t in xrange(len(countPerTemp))])
				p_unb[q1][q2][e] = float(enu) / float(den)
	return p_unb

def punbEffDec(keys, histoPerTemp3Deff, f, ct0, countPerTemp):
	p_unb = {}
	
	for l in keys:  # keys collected from all temperatures
		e = l[2]
		enu = sum([histoPerTemp3Deff[t][l] if l in histoPerTemp3Deff[t] else 0.0  for t in xrange(len(countPerTemp))])
		den = sum([Decimal(countPerTemp[t]) * Decimal(f[t]) * ct0[e][t] for t in xrange(len(countPerTemp))])
		p_unb[l] = float(Decimal(enu) / Decimal(den))
	
	return p_unb

def punbEff(keys, histoPerTemp3Deff, f, c, countPerTemp):
	p_unb = {}

	for l in keys:  # keys collected from all temperatures
		e = l[-1] # energy bin index
		assert e < len(c)
		
		enu = sum([histoPerTemp3Deff[t][l] if l in histoPerTemp3Deff[t] else 0.0  for t in xrange(len(countPerTemp))])

		den = sum([countPerTemp[t] * f[t] * c[e][t] for t in xrange(len(countPerTemp))])
		
		p_unb[l] = float(enu) / float(den)
	
	return p_unb

def KtoPRFbeta(temp):
	return 1.0 / (temp * kelvin_in_pru)

def countPerQuadrantPerT(Tindex, histogram, QAbins, QBbins, qa_bisect, qb_bisect, verbose=True):
	return countPerQuadrant(histogram, [Tindex], QAbins, QBbins, qa_bisect, qb_bisect, verbose)

def countPerQuadrant(histogram, temps, QAbins, QBbins, qa_bisect, qb_bisect, verbose=True):
	if verbose: print "Counts per quadrant:"
	if verbose: print "qa,qb bisects:", qa_bisect, qb_bisect
	if verbose: print "lower left, upper right, upper left, lower right"
	quadrants = [0, 0, 0, 0] 
	quadrantsByTemp = [[0, 0, 0, 0] for t in xrange(len(temps))]

	for t in xrange(len(temps)):
		for qa in xrange(len(QAbins)):
			QAtmp = QAbins[qa]
			for qb in xrange(len(QBbins)):
				count = histogram[t][qa][qb]
			
				QBtmp = QBbins[qb]
				if QAtmp < qa_bisect and QBtmp < qb_bisect:
					quadrants[0] += count
					quadrantsByTemp[t][0] += count
				elif  QAtmp >= qa_bisect and QBtmp >= qb_bisect:
					quadrants[1] += count
					quadrantsByTemp[t][1] += count
				elif  QAtmp >= qa_bisect and QBtmp < qb_bisect:
					quadrants[2] += count
					quadrantsByTemp[t][2] += count
				elif  QAtmp < qa_bisect and QBtmp >= qb_bisect:
					quadrants[3] += count
					quadrantsByTemp[t][3] += count
				else:
					print "error"
	if verbose: print quadrants, "total"
	
	for t in xrange(len(temps)):
		if verbose: print quadrantsByTemp[t], "T", t
	
	return quadrantsByTemp

def Diffusivity(T, f, temps):
	tmp = []
	for i in xrange(len(temps)):
		tmp.append(temps[i][1])
	tmp = list(reversed(tmp))
	for i in xrange(len(tmp) - 1):
		t_lo = tmp[i]
		t_hi = tmp[i + 1]
		
		if T >= t_lo and T < t_hi:
			f_lo = f[i]
			f_hi = f[i + 1]
			dT = t_hi - t_lo
			df = f_hi - f_lo
			dfoverdT = df / dT
			print T, dT, df, dfoverdT
			return dT / dfoverdT
	return None

def popt(T, f, tkel, C):
	for i in xrange(len(tkel) - 1):
		t_lo = tkel[i]
		t_hi = tkel[i + 1]
	
		if T >= t_lo and T < t_hi:
			f_lo = f[i]
			f_hi = f[i + 1]
			dT = t_hi - t_lo
			df = f_hi - f_lo
			
			return C * math.sqrt(abs((1.0 / dT) * (df / dT)))
	sys.exit()

def psumT(T, f, tkel, C):
	s = 0.0
	
	for i in xrange(len(tkel) - 1):
		t_lo = tkel[i]
		t_hi = tkel[i + 1]
		f_lo = f[i]
		f_hi = f[i + 1]
		dT = t_hi - t_lo
		df = f_hi - f_lo
		s += C * math.sqrt(abs((1.0 / dT) * (df / dT))) * dT
		
		if T >= t_lo and T < t_hi:
			break
	return s

def geomTemp(Tstart, Tend, ntmps):
	T = [Tstart]
	for i in xrange(ntmps - 2):
		Tnext = T[i] * math.pow(float(Tend) / Tstart, (1.0 / (ntmps - 1)))
		T.append(Tnext)
	T.append(Tend)
	assert len(T) == ntmps
	return T

def swapProb(temps):
	p = []
	for t in xrange(len(temps)):
		for i in xrange(len(temps)):
			f = math.exp((temps[t] - temps[i]) / temps[t])
			if t < len(temps) - 1 and temps[i] == temps[t + 1]:
				p.append(f)
	return p

def optTbyOscillator(temp, cv, rt=5.0, maxit=5000, adjustRange=None):
	if adjustRange != None:
		newC = []
		for i in xrange(len(adjustRange)):
			T = adjustRange[i]
			nearestTdist = 100000
			nTd_index = None
			for j in xrange(len(temp)):
				if abs(temp[j] - T) < nearestTdist:
					nearestTdist = abs(temp[j] - T)
					nTd_index = j
			
			newC.append(cv[nTd_index])
		
		t = adjustRange
		c = newC
	else:	
		t = temp
		c = cv
	
	assert len(t) == len(c)
	
	c_peaks = []
	ends = [0, len(c) - 1]
	
	for i in range(1, len(c) - 1):
		if c[i - 1] < c[i] and c[i + 1] < c[i]:
			if not i in c_peaks:
				c_peaks.append(i)
	
	dT = []
	for i in xrange(len(t) - 1):
		dT.append(t[i + 1] - t[i])
	
	prevX = [i for i in t]
	
	prevE = 1000000
	for IT in xrange(maxit):
		
		x = [i for i in prevX]
		# perturb
		for d in xrange(len(t)):
			r = uniform(-rt, rt)
			if not d in c_peaks and not d in ends:
				
				x[d] = x[d] + r 
			else:
				x[d] = t[d]
		assert len(x) == len(t)
		dX = []
		for i in xrange(len(x) - 1):
				dX.append(x[i + 1] - x[i])
		
		# score
		
		E = 0.0
		for d in xrange(len(t) - 1):
			dC = c[d + 1] - c[d] / t[d + 1] - t[d]
			aC = (c[d + 1] + c[d]) / 2.0
			
			p1 = max(c) / 1.0  # aC
			p2 = 1.3
			p3 = 2.0
			
			k = p2 * math.pow((abs(dC)) / (p1 + abs(dC)), (1.0 / p3))
			
			if d in c_peaks:
				# left point is peak
				E += math.pow(dX[d] - dT[d] * (1 - k), 2)  # (0.2)
			elif d + 1 in c_peaks:
				# right point is peak
				E += math.pow(dX[d] - dT[d] * (1 - k), 2)  # (0.2)
			elif d - 1 in c_peaks:
				E += math.pow(dX[d] - dT[d] * (1 - k), 2)  # (0.4)
			elif d + 2 in c_peaks:
				E += math.pow(dX[d] - dT[d] * (1 - k), 2)  # (0.4)
			elif d - 2 in c_peaks:
				E += math.pow(dX[d] - dT[d] * (1 - k), 2)  # (0.4)
			elif d + 3 in c_peaks:
				E += math.pow(dX[d] - dT[d] * (1 - k), 2)  # (0.4)
			else:
				# not peak-adjacent
				E += 0.4 * math.pow(dX[d] - dT[d] * (1 - k), 2)  # 0.01
		
		if E < prevE:
			accept = 1.0
		else:
			accept = math.exp(prevE - E)
		
		if accept == 1.0:
			prevE = E
			prevX = [i for i in x]
			
		else:
			r = uniform(0.0, 1.0)
			if accept > r:
				prevE = E
				prevX = [i for i in x]
		
		# print IT,prevE,accept,prevX
	
	# print c_peaks
	return prevX, c_peaks

def interpolate(X, Y):
	f = interp1d(X, Y, kind="cubic")
	return f

def distanceOnCurve(x1, x2, curve, nintervals=10):
	d = 0.0
	interv = (x2 - x1) / float(nintervals)
	x = [x1 + (i * interv) for i in xrange(nintervals)]
	
	for i in xrange(len(x) - 1):
		
		d += math.sqrt(math.pow(x[i + 1] - x[i], 2) + math.pow(curve(x[i + 1]) - curve(x[i]), 2))
	return d

def optTbyTM(Tmin, Tmax, Tm, n):
	halftemp = [math.exp for i in xrange(n / 2)]

def optTbyCurve(temp, cv, rt=1.0, maxit=500, adjustRange=None, polyfit=0):
	if adjustRange != None:
		curve = interpolate(temp, cv)
		newC = []
		for i in xrange(len(adjustRange)):
			T = adjustRange[i]
			
			if T < temp[0]:
				newC.append(cv[0])
			elif T > temp[-1]:
				newC.append(cv[-1])
			else:
				newC.append(curve(T))
		
		t = adjustRange
		c = newC
	else:	
		t = temp
		c = cv
	
	assert len(t) == len(c)
	
	Tm = t[c.index(max(c))]
	
	if polyfit > 0:
		fit = numpy.polyfit(t, c, polyfit)
		tmp = numpy.polyval(fit, t)
		curve = interpolate(t, tmp)
	else:
		curve = interpolate(t, c)
	
	c_peaks = []
	ends = [0, len(c) - 1]
	
	for i in range(1, len(c) - 1):
		if c[i - 1] < c[i] and c[i + 1] < c[i]:
			if not i in c_peaks:
				c_peaks.append(i)
	
	dL = []
	for i in xrange(len(t) - 1):
		
		dL.append(distanceOnCurve(t[i], t[i + 1], curve))  # math.sqrt( math.pow(t[i+1] - t[i],2) + math.pow(c[i+1] - c[i],2) ) )
	print c
	print curve(numpy.array(t))
	totalL = float(sum(dL))
	targetInterval = totalL / (len(dL))
	
	print totalL, targetInterval, len(dL)
	
	prevX = [t[0]]
	newC = [c[0]]
	thresh = 1
	for i in range(1, len(t) - 1):
		go = True
		tmp = prevX[i - 1]
		print tmp
		diff = distanceOnCurve(prevX[i - 1], tmp, curve) 
		while abs(diff - targetInterval) > thresh:
				
			if tmp + 0.01 > t[-1]:
				break
			ctmp = curve(tmp)
			diff = distanceOnCurve(prevX[i - 1], tmp, curve) 
			if diff > targetInterval:
				tmp -= 0.01
			elif diff < targetInterval:
				tmp += 0.01
			else:
				break
			
		print tmp, ctmp, diff, targetInterval
		newC.append(float(ctmp))
		prevX.append(tmp)
	prevX.append(t[-1])
	newC.append(c[-1])
	print prevX, len(prevX)
	print newC, len(newC)
	return prevX, c_peaks, newC

def findMinimumInHisto(histo, QbinsA, QbinsB, minqa, maxqa, minqb, maxqb):
	minimum = 10000000
	minBinA = None
	minBinB = None
	
	binsizeA = QbinsA[1] - QbinsA[0]
	binsizeB = QbinsB[1] - QbinsB[0]
	
	for i in xrange(len(QbinsA)):
		if QbinsA[i] >= minqa and QbinsA[i] + binsizeA <= maxqa:
			for j in xrange(len(QbinsB)):
				if QbinsB[j] >= minqb and QbinsB[j] + binsizeB <= maxqb:
					if histo[i][j] < minimum:
						minimum = histo[i][j]
						minBinA = i
						minBinB = j
	
	return minimum, minBinA, minBinB	

def findMinimumInHisto1D(histo, bins, minq, maxq):
	minimum = 10000000
	minBin = None
	# minBinB = None
	
	binsize = bins[1] - bins[0]
	# binsizeB = QbinsB[1] - QbinsB[0]
	
	for i in xrange(len(bins)):
		if bins[i] >= minq and bins[i] + binsize <= maxq:
				if histo[i] < minimum:
					minimum = histo[i]
					minBin = i
	
	return minimum, minBin	

def tempsFromDF(df):
	temps = {}
	for ti in df["temperature"]:
		if not ti in temps:
			temps[ti] = 1
		else:
			temps[ti] += 1
	sorttemps = sorted(temps.keys())
	tfreqs = [temps[s] for s in sorttemps]
	return sorttemps, tfreqs

def wham(histoPerTemp3Deff, T, countPerTemp, Emidpoints, max_wham_iterations, histokeys, temps, useDecimalPackage=False, verbose=False, wham_cut=1e-3):
	nEtotbins = len(Emidpoints)
	if useDecimalPackage:
		Decimal.getcontext().prec = 9
		c = [[     Decimal(0.0)     for t in xrange(len(temps))] for e in xrange(nEtotbins)]  
		
		for e in xrange(nEtotbins):
			for t in xrange(len(temps)):
				c[e][t] = Decimal(-(temps[t] - (1.0 / T)) * Emidpoints[e]).exp() 
	else:
		c = [[     math.exp(-(temps[t] - (1.0 / T)) * Emidpoints[e])     for t in xrange(len(temps))] for e in xrange(nEtotbins)]  

	p_unb = 0

	p_old = copy.copy(p_unb)
	done = False
	f = [1.0 for t in xrange(len(temps))]
	iterations = 0
	while(not done):
		iterations += 1
		if useDecimalPackage:
			p_unb = punbEffDec(histokeys, histoPerTemp3Deff, f, c, countPerTemp)
			f = [ float(Decimal(1.0) / sum([ c[l[2]][t] * Decimal(p_unb[l]) for l in p_unb ]))  for t in xrange(len(temps))]
		else:
			p_unb = punbEff(histokeys, histoPerTemp3Deff, f, c, countPerTemp)
			f = [ float(1.0 / sum([ c[l[-1]][t] * p_unb[l] for l in p_unb ]))  for t in xrange(len(temps))]
		
		if  not max_wham_iterations == -1 and iterations > max_wham_iterations: 
			done = True
		
		
		if p_old != None and p_unb != None:
			p_delta = max([ abs(p_old[l] - p_unb[l]) for l in histokeys if l in p_unb and l in p_old])
		
			print p_delta, wham_cut
		
			if p_delta <= wham_cut:
				done = True
		
		p_old = copy.copy(p_unb)
			
	return p_unb

def polyfitCv(tprf, cvs, polyfitDegree=10, maxTmCut=2, minTmCut=0):
	# test best fit
	deviations = {}
	mindev = 999e999
	mindev_deg = -1
	for deg in range(1, polyfitDegree + 1):
		fit = np.polyfit(tprf, cvs, deg)
		
		testcv = np.polyval(fit, tprf)
		rmsd = math.sqrt(sum([math.pow(cvs[i] - testcv[i], 2) for i in xrange(len(cvs))]))
		
		deviations[deg] = rmsd
		if rmsd < mindev: 
			mindev = rmsd
			mindev_deg = deg
	
	print mindev, mindev_deg
	fit = np.polyfit(tprf, cvs, mindev_deg)
	
	finegrid = np.arange(min(tprf), max(tprf), (max(tprf) - min(tprf)) / 100.0) 
	fineCv = np.polyval(fit, finegrid)
	
	dCv = [ fineCv[i + 1] - fineCv[i - 1]  if i > 0 and i < len(fineCv) - 1 else 0 for i in xrange(len(fineCv))]
	dCv[0] = dCv[1]
	dCv[-1] = dCv[-2]
	dCv = [ dCv[i + 1] - dCv[i - 1]  if i > 0 and i < len(dCv) - 1 else 0 for i in xrange(len(dCv))]
	dCv[0] = dCv[1]
	dCv[-1] = dCv[-2]
	maxcv = max([fineCv[c] for c in xrange(len(fineCv)) if finegrid[c] >= minTmCut and finegrid[c] <= maxTmCut])
	maxcv_i = list(fineCv).index(maxcv)
	Tmelti = maxcv_i
	Tmelt = finegrid[Tmelti]
	cvs = list(fineCv)
	return finegrid, cvs, Tmelt, maxcv

def binning(col, break_points, labels=None):
	if not labels:
		labels = range(len(break_points)-1)
	#print labels
	#Binning using cut function of pandas
	#print break_points, len(break_points)
	colBin = pd.cut(col, bins=break_points, labels=labels, precision=6, include_lowest=True)
	#print "colBin:",colBin
	#print colBin.shape
	#print "bins:",bins, len(bins)
	#sys.exit()
	assert not colBin.hasnans()
	return colBin

def get3Dhisto(df, ai, bi, QAbins, QBbins, Etotbins, temps, mintime=0, maxtime=-1,ci="Time_in_MC_cycles", ti="temperature", ei="Etot"):
	histo3DPerTemp = [{} for t in xrange(len(temps))]
	histokeys = []
	if maxtime==-1:
		maxtime = df[ci].max()
	
	
	A = df[ai]
	B = df[bi]
	E = df[ei]

	
	Abreaks = QAbins+[2*QAbins[-1] - QAbins[-2]]
	Bbreaks = QBbins+[2*QBbins[-1] - QBbins[-2]]
	Ebreaks = Etotbins+[2*Etotbins[-1] - Etotbins[-2]]
	
	
	assert not A.hasnans()
	assert Abreaks[0] == min(Abreaks) and Abreaks[-1] == max(Abreaks)
	assert all([a >= Abreaks[0]-1e-10 and a<=Abreaks[-1]+1e-10 for a in A ])
	
	assert not B.hasnans()
	assert Bbreaks[0] == min(Bbreaks) and Bbreaks[-1] == max(Bbreaks)
	assert all([b >= Bbreaks[0]-1e-10 and b <= Bbreaks[-1]+1e-10 for b in B ])
	
	assert not E.hasnans()
	assert Ebreaks[0] == min(Ebreaks) and Ebreaks[-1] == max(Ebreaks)
	assert all([e >= Ebreaks[0]-1e-10 and e<=Ebreaks[-1]+1e-10 for e in E ])
	
	assert len(A)==len(B)==len(E)
	
	for t in xrange(len(temps)):
		tbin = temps[t]
		tempmatch = df[ti] == tbin
		abovemin = df[ci] >= mintime
		belowmax = df[ci] <= maxtime

		AtempBin = binning(A[tempmatch & abovemin & belowmax], Abreaks)
		BtempBin = binning(B[tempmatch & abovemin & belowmax], Bbreaks)
		EtempBin = binning(E[tempmatch & abovemin & belowmax], Ebreaks)
		
		keys = zip(AtempBin, BtempBin, EtempBin)
		
		for k in keys:
			assert all(type(i) in [int, np.int64, np.int32] for i in k), repr(k)
			assert k[0]>=0 and k[0]<=len(Abreaks)-1, k
			assert k[1]>=0 and k[1]<=len(Bbreaks)-1, k
			assert k[2]>=0 and k[2]<=len(Ebreaks)-1, k
			if k in histo3DPerTemp[t]:
				histo3DPerTemp[t][k] += 1
			else:
				histo3DPerTemp[t][k] = 1
				histokeys.append(k)
		
	return histo3DPerTemp, histokeys

def specificHeat(temps, etotPerTemp, maxTmCut=2, minTmCut=0):
	
	cvs = [ (1.0 / math.pow(t, 2)) * (np.mean([math.pow(i, 2) for i in etotPerTemp[t]]) - math.pow(np.mean([i for i in etotPerTemp[t]]), 2)) for t in temps ]
	
	tmp = [cvs[c] for c in xrange(len(cvs)) if temps[c] <= maxTmCut and temps[c] >= minTmCut ]
	
	#dCv = [cvs[i] - cvs[i - 1] if i > 0 else 0 for i in xrange(len(cvs))]
	
	
	Tmelt = None
	maxcv = max(tmp)
	for i in xrange(len(cvs)):
		if cvs[i] == maxcv:
			Tmelt = temps[i]
	
	return cvs, Tmelt, maxcv

def getHistograms(raw_data, QAbins, QBbins, Etotbins, temps, verbose=False):
	histoPerTempQAE = [ {} for t in xrange(len(temps))]
	histoPerTempQBE = [ {} for t in xrange(len(temps))]
	histoPerTempEtot = [ [0 for i in xrange(len(Etotbins))] for t in xrange(len(temps))]
	histoPerTempQ1 = [ [0 for i in xrange(len(QAbins))] for t in xrange(len(temps))]
	histoPerTempQ2 = [ [0 for i in xrange(len(QBbins))] for t in xrange(len(temps))]
	countPerTemp = [0 for t in xrange(len(temps))]
	histoPerTemp3Deff = [{} for t in xrange(len(temps))]
	totalCount = 0
	histokeys = []
	histokeysA = []
	histokeysB = []
	
	for row in raw_data:
		Ti = row[1]
		Etot = row[2]
		QA = row[3]
		QB = row[4]
	
		QAi = getBinIndex(QAbins, QA)
		QBi = getBinIndex(QBbins, QB)
		Etoti = getBinIndex(Etotbins, Etot)
		
		histoPerTempEtot[Ti][Etoti] += 1
		histoPerTempQ1[Ti][QAi] += 1
		histoPerTempQ2[Ti][QBi] += 1
		countPerTemp[Ti] += 1
		
		key = (QAi, QBi, Etoti)
		keyA = (QAi, Etoti)
		keyB = (QBi, Etoti)
		
		histokeys.append(key)
		histokeysA.append(keyA)
		histokeysB.append(keyB)
		
		if keyA in histoPerTempQAE[Ti]:
			histoPerTempQAE[Ti][keyA] += 1
		else:
			histoPerTempQAE[Ti][keyA] = 1
		
		if keyB in histoPerTempQBE[Ti]:
			histoPerTempQBE[Ti][keyB] += 1
		else:
			histoPerTempQBE[Ti][keyB] = 1
		
		if key in histoPerTemp3Deff[Ti]:
			histoPerTemp3Deff[Ti][key] += 1
		else:
			histoPerTemp3Deff[Ti][key] = 1
		totalCount += 1
		
	
	if verbose:
		print "total data row count:", totalCount
		print "total counts per T:", [countPerTemp[t] for t in xrange(len(temps))]
		#countPerQuadrant(histoPerTemp2D, temps, QAbins, QBbins, qa_bisect, qb_bisect)
	
	return histokeys, histokeysA, histokeysB, countPerTemp, histoPerTemp3Deff, totalCount, histoPerTempQAE, histoPerTempQBE , histoPerTempEtot, histoPerTempQ1, histoPerTempQ2

def parseDataFromSingleFile(datafile, temps, indexA, indexB, minsteps, maxsteps, stepsize, nativenessQ, EAmin, EBmin, constTempStats, Tmix=None, weightA=1.0, weightGO=1.0, lundVsGo=False):
	Etotmin = 99999999
	Etotmax = -99999999

	maxtime = 0
	
	raw_data = []

	tempCountPerReplica = {}

	ncols = None

	replica = -1
	prevtime = -1
	raw_data_per_replica = {}

	etotPerTemp = [[] for t in xrange(len(temps))]
	egoPerTemp = [[] for t in xrange(len(temps))]
	lvgPerTemp = [[] for t in xrange(len(temps))]
	nT0 = [0 for t in xrange(len(temps))]
	nT1 = [0 for t in xrange(len(temps))]
	current = {}
	timeStep = None
	firsttimestep = None
	
	stepcounter = stepsize
	
	for line in open(datafile).readlines() :
		tmp = line.replace("\x00", "").strip().split()
	
		if ncols == None:
			ncols = len(tmp)
		if len(tmp) != ncols: 
			print "T=", tmp[0], " -> mismatch of column number!", len(tmp), " vs. ", ncols
			print tmp
			continue
		if tmp[0] == '': print tmp
		try:
			time = int(tmp[0])
		except ValueError:
			print tmp
			continue
		
		if firsttimestep == None:
			firsttimestep = time
		
		if prevtime == -1:
			replica = 0
			current[replica] = -1
		elif time == firsttimestep:
			replica += 1
			current[replica] = -1
		else:
			if timeStep == None:
				timeStep = time - prevtime
		prevtime = time
		
		if minsteps > 0 and time < minsteps:
			continue
		if maxsteps > 0 and time > maxsteps:
			continue
		
		if time > maxtime: maxtime = time
		
		if stepcounter != 1:
		
			stepcounter -= 1
		
			if stepcounter >= 0: 
			
				continue
			else:
				stepcounter = stepsize
			
		Ti = int(tmp[1])
		if constTempStats:
			if Ti != 0:
				print "WARNING: at constant temperature there should only be index 0, but found ", Ti
				Ti = 0
		
		Etot = float(tmp[2])
		
		exvol = float(tmp[3])
		locexvol = float(tmp[4])
		bias = float(tmp[5])
		torsion = float(tmp[6])
		hbmm = float(tmp[7])
		hbms = float(tmp[8])
		hydrop = float(tmp[9])
		charged = float(tmp[10])
		
		Etrans = exvol + locexvol + bias + torsion + hbmm + hbms + hydrop + charged
		
		EA = float(tmp[indexA])
		EB = float(tmp[indexB])
		
		Eatt = hbmm + hbms + hydrop + charged
		
		Ego = Etot - Etrans
		
		newEA = EA
		newEB = EB
		
		newEAmin = EAmin
		newEBmin = EBmin
		
		newEgo = Ego
		
		if weightA != 1.0 or weightGO != 1.0:
			if Tmix == None:
				assert EA + EB - Ego < 0.0001, (EA + EB - Ego)
				newEA = EA * weightA * weightGO
				newEB = EB * weightGO
			
				newEAmin = EAmin * weightA * weightGO
				newEBmin = EBmin * weightGO
			
				newEgo = newEA + newEB
			else:
				newEgo = -Tmix * math.log(math.exp((-EA * weightA * weightGO) / Tmix) + math.exp((-EB * weightGO) / Tmix))
		
		Etot_old = Etot
		Etot = Etrans + newEgo
		
		if Ti == 0:
			if current[replica] == -1: current[replica] = 0
			else:
				if current[replica] == 1:
					current[replica] = 0
		elif Ti == len(temps) - 1:
			if current[replica] == -1: current[replica] = 1
			else:
				if current[replica] == 0:
					nT1[Ti] += 1
					current[replica] = 1
	
		if current[replica] == 0:
			nT0[Ti] += 1
		elif current[replica] == 1:
			nT1[Ti] += 1
	
		etotPerTemp[Ti].append(Etot)
		lvgPerTemp[Ti].append(newEgo - Eatt)
		egoPerTemp[Ti].append(newEgo)
	
		if replica in tempCountPerReplica:
			if Ti in tempCountPerReplica[replica]:
				tempCountPerReplica[replica][Ti] += 1
			else:
				tempCountPerReplica[replica][Ti] = 1
		else:
			tempCountPerReplica[replica] = {}
	
		if Etot < Etotmin: Etotmin = Etot
		if Etot > Etotmax: Etotmax = Etot
	
		if nativenessQ in ["True", "true", "1"]:
			QA = float(tmp[indexA])
			QB = float(tmp[indexB])
		else:	
			QA = -float(newEA) / newEAmin
			QB = -float(newEB) / newEBmin
			if QA < 0.0: QA = 0.0
			if QB < 0.0: QB = 0.0
		
		if replica in raw_data_per_replica:
			raw_data_per_replica[replica].append([time, Ti, Etot, QA, QB])
		else:
			raw_data_per_replica[replica] = [ [time, Ti, Etot, QA, QB] ]
		
		raw_data.append([time, Ti, Etot, QA, QB])
	
	if maxtime < minsteps:
		print "ERROR: no data - lower minsteps!"
		sys.exit(1)
	
	if lundVsGo:
		for t in xrange(len(temps)):
			print temps[t][1], numpy.mean(lvgPerTemp[t]), numpy.std(lvgPerTemp[t]), min(lvgPerTemp[t]), max(lvgPerTemp[t])
	
	return raw_data, etotPerTemp, egoPerTemp, nT0, nT1, tempCountPerReplica, Etotmin, Etotmax, maxtime, raw_data_per_replica, timeStep

def getCv(temps, etotPerTemp, maxTmCut=2000, minTmCut=0, polyfitDegree=0, Tunit="Profasi"):
	if Tunit == "Kelvin":
		tkel = [temps[t][1] for t in reversed(xrange(len(temps)))]
	else:
		tkel = [temps[t][0] for t in reversed(xrange(len(temps)))]
		tprf = [temps[t][0] for t in reversed(xrange(len(temps)))]
	
	tprf = [temps[t][0] for t in reversed(xrange(len(temps))) ]
	print temps
	cvs = [ (1.0 / math.pow(tprf[t], 2)) * (numpy.mean([math.pow(i, 2) for i in etotPerTemp[t]]) - math.pow(numpy.mean([i for i in etotPerTemp[t]]), 2))   for t in reversed(xrange(len(temps))) ]
	
	tmp = [cvs[c] for c in xrange(len(cvs)) if tprf[c] <= maxTmCut and tprf[c] >= minTmCut ]
	
	dCv = [cvs[i] - cvs[i - 1] if i > 0 else 0 for i in xrange(len(cvs))]
	
	if polyfitDegree >0:
		# test best fit
		deviations = {}
		mindev = 999e999
		mindev_deg = -1
		for deg in range(1, polyfitDegree + 1):
			fit = numpy.polyfit(tprf, cvs, deg)
			
			testcv = numpy.polyval(fit, tprf)
			rmsd = math.sqrt(sum([math.pow(cvs[i] - testcv[i], 2) for i in xrange(len(cvs))]))
			
			deviations[deg] = rmsd
			if rmsd < mindev: 
				mindev = rmsd
				mindev_deg = deg
		
		print mindev, mindev_deg
		fit = numpy.polyfit(tprf, cvs, mindev_deg)
		
		finegrid = numpy.arange(min(tprf), max(tprf), (max(tprf) - min(tprf)) / 100.0) 
		fineCv = numpy.polyval(fit, finegrid)
		
		dCv = [ fineCv[i + 1] - fineCv[i - 1]  if i > 0 and i < len(fineCv) - 1 else 0 for i in xrange(len(fineCv))]
		dCv[0] = dCv[1]
		dCv[-1] = dCv[-2]
		dCv = [ dCv[i + 1] - dCv[i - 1]  if i > 0 and i < len(dCv) - 1 else 0 for i in xrange(len(dCv))]
		dCv[0] = dCv[1]
		dCv[-1] = dCv[-2]
		maxcv = max([fineCv[c] for c in xrange(len(fineCv)) if finegrid[c] >= minTmCut and finegrid[c] <= maxTmCut])
		maxcv_i = list(fineCv).index(maxcv)
		Tmelti = maxcv_i
		Tmelt = finegrid[Tmelti]
		cvs = list(fineCv)
	else:
		Tmelt = None
		Tmelti = None
		maxcv = max(tmp)
		for i in xrange(len(cvs)):
			if cvs[i] == maxcv:
				Tmelt = tprf[i]
				Tmelti = i
	
	return cvs, Tmelt, Tmelti, maxcv, dCv

def getDiffusivity(nT0, nT1, temps):
	tkel = [temps[t][1] for t in reversed(xrange(len(temps)))]
	
	fT = list(reversed([ nT1[i] / float(nT0[i] + nT1[i]) if float(nT0[i] + nT1[i]) != 0.0 else 0.0 for i in xrange(len(temps))]))
	
	x = 1.0

	prob = [x * math.sqrt(abs((1.0 / (tkel[i + 1] - tkel[i])) * (fT[i + 1] - fT[i]) / (tkel[i + 1] - tkel[i]))) for i in range(0, len(tkel) - 1)]
	psum = sum([ prob[i] * (tkel[i + 1] - tkel[i]) for i in range(0, len(tkel) - 1)])
	# print psum
	
	if psum != 0.0:
		x = 1.0 / psum
	else:
		x = 0.0
	prob = [x * math.sqrt(abs((1.0 / (tkel[i + 1] - tkel[i])) * (fT[i + 1] - fT[i]) / (tkel[i + 1] - tkel[i]))) for i in range(0, len(tkel) - 1)]
	psum = sum([ prob[i] * (tkel[i + 1] - tkel[i]) for i in range(0, len(tkel) - 1)])
	Diff = [ math.log(math.pow(1.0 / prob[i], 2)) if prob[i] != 0.0 else None for i in xrange(len(prob))]
	
	return fT, prob, Diff

def expfunc(x, a, b=1.0, c=0.0):
    return a * numpy.exp(-b * x) + c

def linfunc(x, a, b):
	return a * x + b

def gaussfunc(x, a, b):
	return a * numpy.exp(-b * numpy.power(x, 2.0))

def gaussfuncScalar(x, a, b, m=1.0):
	return a * math.exp(-b * math.pow(x - m, 2.0))

def autoCorrelation(temps, etotPerTemp, datafile, maxtime, tau_interval=1, plots=True, nSweepsPerData=1000, acceptACcut=0.01, trajfracshown=8, showPlots=False, filename="ac"):
	print "temps=%i, etotPerTemp=[%i][%i], tau_interval=%i, plots=%i, nSweepsPerData=%i, acceptACcut=%f, trajfracshown=%i" % (len(temps), len(etotPerTemp), len(etotPerTemp[0]), tau_interval, plots, nSweepsPerData, acceptACcut, trajfracshown)
	ac = {}
	popts = []
	tauAccept = [None for t in xrange(len(temps))]
	for t in xrange(len(temps)):
		E = etotPerTemp[t]
		ac[t] = []
		meanE = numpy.mean(E)
		stdE = numpy.std(E)
		taus = range(1, len(E) / trajfracshown, tau_interval)
		for tau in taus:
			eterm = numpy.mean([ (E[i + tau] - meanE) * (E[i] - meanE)   for i in xrange(len(E) - tau) ]) / math.pow(stdE, 2)
			
			ac[t].append(eterm)
			if tauAccept[t] == None and eterm < acceptACcut:
				tauAccept[t] = tau
				print "T=%s - AC(tau=%i) < %s; " % (str(round(temps[t][0], 3)), tau, str(acceptACcut)),
		
		popt = None
		try:
			assert len(taus) == len(ac[t])
			tmp_taus = numpy.array([taus[ti] for ti in xrange(len(taus)) if taus[ti] <= 100])
			tmp_ac = numpy.array([ac[t][ti] for ti in xrange(len(taus)) if taus[ti] <= 100])
			linpopt, linpcov = curve_fit(linfunc, tmp_taus , tmp_ac)
			popt = linpopt
			print " fit: f(x)=ax+b; [a,b]=[%f,%f]; -b/a=%f" % (popt[0], popt[1], -popt[1] / popt[0])
			popts.append(-popt[1] / popt[0])
		except RuntimeError as re:
			print re
			popts.append(None)
			
		color = 0.0 + (float(t) / (len(temps) + 1))
		
		if plots: 
			p.plot(taus, ac[t], label=str(round(temps[t][0], 2)), linewidth=1, color=str(color))
			p.scatter(popts , [0 for i in xrange(len(popts))], marker="x", color="0.5")
	if plots: 
		if temps[0] != temps[1]: p.legend(ncol=3)
		p.title(os.path.basename(datafile).replace("_rt.txt", "") + ";T%i;N%i" % (maxtime, len(etotPerTemp[0])))
		# p.yscale('log')
		p.xlabel("tau (=%i MC sweeps)" % nSweepsPerData)
		p.ylabel("auto correlation (tau)")
		p.grid(b=True, which='both')
		p.savefig(filename + "_tauVSac.png")
		if showPlots: p.show()
	kelvins = [temps[t][0] for t in xrange(len(temps))]
	
	means = [ numpy.mean([j   for j in ac[t] if j != None]) for t in xrange(len(temps))]
	print means
	
	cvs, Tmelt, Tmelti, maxcv, dCv = getCv(temps, etotPerTemp)
	if plots: 
		fig, ax1 = p.subplots()
		ax1.set_title(os.path.basename(datafile).replace("_rt.txt", "") + ";T%i;N%i" % (maxtime, len(etotPerTemp[0])))
		ax1.set_xlabel("T")
		ax1.set_ylabel("tau | AC(tau) < %s" % str(acceptACcut))
		ax1.plot(kelvins, tauAccept, marker="x")
		ax2 = ax1.twinx()
		ax2.plot(kelvins, cvs, marker="", color="k", linewidth=1, linestyle=":")
		ax2.set_ylabel("Cv(T)")
		p.savefig(filename + "_tauVcv.png")
		if showPlots: p.show()
	
	return list(reversed(means))

def doConstTempStats(raw_data_per_replica):
	maxt = [0 for r in xrange(len(raw_data_per_replica))]
	maxqa = [0 for r in xrange(len(raw_data_per_replica))]
	maxqb = [0 for r in xrange(len(raw_data_per_replica))]
	
	for ri in xrange(len(raw_data_per_replica)):
		if not ri in raw_data_per_replica:
			print "missing replica:", ri
			continue
		r = raw_data_per_replica[ri]
		for rj in xrange(len(r)):
			time = int(r[rj][0])
			qa = float(r[rj][3])
			qb = float(r[rj][4])
		
			if time > maxt[ri]: 
				maxt[ri] = time
				maxqa[ri] = qa
				maxqb[ri] = qb
	
	counter = [0, 0, 0, 0, 0, 0]
	for ri in xrange(len(raw_data_per_replica)):
		if maxqa[ri] <= 0.6 and maxqb[ri] <= 0.6:
			 counter[0] += 1
		elif maxqa[ri] > 0.6 and maxqb[ri] > 0.6:
			 counter[1] += 1
		elif maxqa[ri] > 0.6 and maxqb[ri] <= 0.6:
			counter[2] += 1
		elif maxqa[ri] <= 0.6 and maxqb[ri] > 0.6:
			counter[3] += 1
		if maxqa[ri] > 0.6: 
			 counter[4] += 1
		if maxqb[ri] > 0.6:
			counter[5] += 1
	if not float(sum(counter[:4])) == 0:
		
		print "QA <= 0.6 and QB <= 0.6 : %s%%" % (100 * counter[0] / float(sum(counter[:4])))
		print "QA >  0.6 and QB >  0.6 : %s%%" % (100 * counter[1] / float(sum(counter[:4])))
		print "QA >  0.6 and QB <= 0.6 : %s%%" % (100 * counter[2] / float(sum(counter[:4])))
		print "QA <= 0.6 and QB >  0.6 : %s%%" % (100 * counter[3] / float(sum(counter[:4])))
		print "QA >  0.6               : %s%%" % (100 * counter[4] / float(sum(counter[:4])))
		print "QB >  0.6               : %s%%" % (100 * counter[5] / float(sum(counter[:4])))
	else:
		print "No data available!"
	return counter

def sumUnderCurve(x, y):
	s = 0.0
	dxs = []
	dys = []
	
	assert len(x) == len(y)
	
	for i in xrange(len(x) - 1):
		if y[i] != None and y[i + 1] != None:
			dX = x[i + 1] - x[i]
			dxs.append(dX)
			dY = (y[i + 1] + y[i]) / 2.0
			dys.append(dY)
			s += dX * dY 
	
	return s, dxs, dys

def cooperativity(Tkel, Cv, Tm, Cvmax, left_end_index=None, right_begin_index=None):
	print "cooperativity.."
	assert len(Tkel) == len(Cv)
	tmi = Tkel.index(Tm)
	assert tmi == Cv.index(Cvmax)
	
	linCvcut = Cvmax / 20.0
	print "linCvcut:", linCvcut
	
	pre_s, pre_dx, pre_dy = sumUnderCurve(Tkel, Cv)
	
	if left_end_index == None:
		left_end_index = len(Cv) - 1
		# left arm
		for dyi in xrange(len(pre_dy[:tmi]) - 1):
			if pre_dy[dyi] / Cvmax > linCvcut:
				left_end_index = dyi
				break
	
	left_lin_popt, left_lin_pcov = curve_fit(linfunc, numpy.array(Tkel[:left_end_index]) , numpy.array(Cv[:left_end_index]))
	print Tkel[:left_end_index]
	print Cv[:left_end_index]
	print "left_lin_popt", left_lin_popt
	print left_lin_pcov
	if right_begin_index == None:
		right_begin_index = 0
		# right arm
		for dyi_tmp in xrange(len(pre_dy[tmi:]) - 1):
			dyi = dyi_tmp + tmi
			if pre_dy[dyi] / Cvmax < linCvcut:
				right_begin_index = dyi
				break
	
	print "fitting curve.."
	right_lin_popt, right_lin_pcov = curve_fit(linfunc, numpy.array(Tkel[right_begin_index:]) , numpy.array(Cv[right_begin_index:]))
	print Tkel[right_begin_index:] 
	print Cv[right_begin_index:]
	print "right_lin_popt", right_lin_popt
	print right_lin_pcov
	
	newCv = []
	left_linCv = []
	right_linCv = []
	for i in xrange(len(Cv)):
		leftlin = linfunc(Tkel[i], left_lin_popt[0], left_lin_popt[1])
		left_linCv.append(leftlin)
		if i < tmi:
			newCv.append(Cv[i] - leftlin)
	
	for i in xrange(len(Cv)):
		rightlin = linfunc(Tkel[i], right_lin_popt[0], right_lin_popt[1])
		right_linCv.append(rightlin)
		if i >= tmi:
			newCv.append(Cv[i] - rightlin)
	
	assert len(Cv) == len(newCv)
	
	print "sum under curve.."
	post_s, post_dx, post_dy = sumUnderCurve(Tkel[left_end_index:right_begin_index], newCv[left_end_index:right_begin_index])
	
	calorimetric_enthalpy = post_s
	
	vantHoff_enthalpy = 2 * math.sqrt(Tm * Tm * Cvmax)
	print "vantHoff_enthalpy =", vantHoff_enthalpy
	print "calorimetric_enthalpy =", calorimetric_enthalpy
	return vantHoff_enthalpy / calorimetric_enthalpy, newCv, left_linCv, right_linCv

def isInState(q1, q2, s):
	if q1 >= s[0] and q1 < s[1] and q2 >= s[2] and q2 < s[3]:
		return True
	else:
		return False

def getCoordPDB(pdb, pdb2data, indexA, indexB):
	tmp = pdb2data[pdb][1:]
	# print tmp
	if nativenessQ in ["True", "true", "1"]:
		QA = float(tmp[indexA])
		QB = float(tmp[indexB])
		EA = None
		EB = None
	else:	
		EA = float(tmp[indexA])
		EB = float(tmp[indexB])
		QA = -float(EA * weightA) / (EAmin * weightA)
		QB = -float(EB) / EBmin
		if QA < 0.0: QA = 0.0
		if QB < 0.0: QB = 0.0
	return QA, QB

def checkHeader(datafile):
	labels = []
	f = open(datafile)
	first = f.readline()
	if "Time_in_MC_cycles" in first:
		labels = first.strip().split()
	
	mins = [10e10 for i in xrange(len(labels))]
	maxs = [-10e10 for i in xrange(len(labels))]
	avg = [[] for i in xrange(len(labels))]
	
	for line in f.readlines():
		tmp = line.strip().split()
		if len(tmp) < 8: continue
		try:
			tmp = map(float, tmp)
		except ValueError as ve:
			print tmp
			print ve
			continue
		if len(tmp) == len(labels):
			for j in xrange(len(tmp)):
				if tmp[j] < mins[j]: mins[j] = tmp[j]
				if tmp[j] > maxs[j]: maxs[j] = tmp[j]
				avg[j].append(tmp[j])
	
	return labels, mins, maxs, [numpy.mean(avg[l]) for l in xrange(len(labels))], [numpy.std(avg[l]) for l in xrange(len(labels))]
	
def isPointInRect(px, py, x1, x2, y1, y2):
	if px < x2 and px > x1 and py < y2 and py > y1:
		return True
	else:
		return False 
	
def isPointInRectangles(px, py, rects):
	i = -1
	for x1, x2, y1, y2, c, l in rects:
		assert x2 >= x1 and y2 >= y1
		i += 1
		
		if isPointInRect(px, py, x1, x2, y1, y2):
			return i
		
	return -1

def AA3to1(three_letter_aa):
	AA_3_to_1 = {
					'ALA':'A',
					'ARG':'R',
					'ASN':'N',
					'ASP':'D',
					'CYS':'C',
					'GLN':'Q',
					'GLU':'E',
					'GLY':'G',
					'HIS':'H',
					'ILE':'I',
					'LEU':'L',
					'LYS':'K',
					'MET':'M',
					'PHE':'F',
					'PRO':'P',
					'SER':'S',
					'THR':'T',
					'TRP':'W',
					'TYR':'Y',
					'VAL':'V'}
	return AA_3_to_1[three_letter_aa]

def AA1to3(one_letter_aa):
	AA_1_to_3 = {
					'G':'GLY',
					'A':'ALA',
					'V':'VAL',
					'L':'LEU',
					'I':'ILE',
					'C':'CYS',
					'M':'MET',
					'F':'PHE',
					'Y':'TYR',
					'W':'TRP',
					'P':'PRO',
					'S':'SER',
					'T':'THR',
					'N':'ASN',
					'Q':'GLN',
					'D':'ASP',
					'E':'GLU',
					'H':'HIS',
					'K':'LYS',
					'R':'ARG'}
	return AA_1_to_3[one_letter_aa]

def atomsFromPDB(filename):
	atoms = {}
	for line in open(filename).readlines():
		if line[:4] == "ATOM":
			a = parseATOMline(line)
			atoms[a["anum"]] = a
	return atoms

def atomsFromGRO(atoms_gro_fname):
	atoms_gro = {}
	gro = open(atoms_gro_fname)
	title = gro.readline()
	natoms = int(gro.readline())
	for line in gro.readlines():
		
		if len(line) >= 44:

			resnr = line[:5].strip()
			residue = line[5:10].strip()
			atype = line[10:15].strip()
			anr = line[15:20].strip()
			x = line[20:28].strip()
			y = line[28:36].strip()
			z = line[36:44].strip()
			atoms_gro[int(anr)] = {"resnr":int(resnr), "residue":residue.strip(), "atype":atype.strip(), "x":float(x), "y":float(y), "z":float(z)}
		
	assert natoms == len(atoms_gro.keys())
	return atoms_gro

def translateAtoms(a1, a2, trans_atoms, ref_atoms):
	trans = []
	# print a1,a2
	for a in [a1, a2]:
		trans_aname = trans_atoms[a]["aname"]
		trans_rname = trans_atoms[a]["rname"]
		trans_rnum = trans_atoms[a]["rnum"]
		trans_chain = trans_atoms[a]["chain"]
		found = False
		for ra in ref_atoms.keys():
			
			ref_aname = ref_atoms[ra]["aname"]
			ref_rname = ref_atoms[ra]["rname"]
			ref_rnum = ref_atoms[ra]["rnum"]
			ref_chain = ref_atoms[ra]["chain"]
			
			if ref_aname == trans_aname and ref_rnum == trans_rnum and ref_chain == trans_chain:
				
				found = True
				# print "found atom",a
				break
		if found:
			trans.append(ra)
		else:
			trans.append(None)
	
	return trans[0], trans[1]

def filterSuperEnsembleContacts(models_contact_lists, models_atoms, atypes, ensembleCutoff, contactType, verbose=False):
	if verbose: print "super ensemble contacts"
	contactsDict = {}
	
	atomtypes = copy.copy(atypes)
	
	if contactType == "HYBRID" and not "CA" in atomtypes:
		atomtypes += ['CA']
	
	nModels = len(models_contact_lists)
	contactsPerModel = [[] for i in xrange(nModels)]
	
	for m in xrange(nModels):
		model_contacts = models_contact_lists[m]
		if verbose: print "%i contacts in model %s" % (len(model_contacts), m)
		
		modelatoms = models_atoms[m]
		if verbose: print len(modelatoms), " atoms in modelatoms"
		
		for c in model_contacts:
			a1 = c[0]
			a2 = c[1]
			d = c[2]
			assert modelatoms[a1]["aname"] in atypes, modelatoms[a1]
			assert modelatoms[a2]["aname"] in atypes, modelatoms[a2]
			atype1 = modelatoms[a1]["aname"]
			atype2 = modelatoms[a2]["aname"]
			rnum1 = modelatoms[a1]["rnum"]
			rnum2 = modelatoms[a2]["rnum"]
			if (rnum1, rnum2, atype1, atype2) in contactsDict:
				contactsDict[(rnum1, rnum2, atype1, atype2)].append(d)
			else:
				contactsDict[(rnum1, rnum2, atype1, atype2)] = [d]
	
	# which contacts are found in all models? what are the means and standard deviations of distances?

	counter = 0
	
	finalContacts = []
	dists = []
	
	global_min_dist = 10000
	global_max_dist = -1
	min_range = 100000
	max_range = -1
	
	if verbose: print "Collected %i contacts from all models." % len(contactsDict.keys())
	
	for i in sorted(contactsDict.keys()):
		
		contactCount = len(contactsDict[i]) 
		assert contactCount <= nModels, str([contactCount, nModels])
		
		frac = contactCount / float(nModels)
		
		if frac >= ensembleCutoff:  # in how many models was this contact found ? Is the fraction above threshold?
			counter += 1
			
			min_d = min(contactsDict[i])
			max_d = max(contactsDict[i])
			range_d = max_d - min_d
			
			if min_d < global_min_dist: global_min_dist = min_d
			if max_d > global_max_dist: global_max_dist = max_d
			if range_d < min_range: min_range = range_d
			if range_d > max_range: max_range = range_d
			
			if verbose: print "(%s) accepted atom pair:" % str(round(frac, 1)), i, " - ", contactCount, "models with this contact.", "Min,Max:", min_d, max_d, " Range:", range_d
			
			if i[2] in atomtypes and i[3] in atomtypes:
				fc = list(i) + [contactsDict[i]]
				finalContacts.append(fc)
			else:
				if verbose: print "atom is not allowed for contact!", ref_atoms[i[0]]["aname"], ref_atoms[i[1]]["aname"]
		else:
			if verbose: print "(%s) REJECTED atom pair:" % str(round(frac, 1)), i, " - ", contactCount, "models with this contact."
	
	if verbose: 
		print "Min and Max distances within all contacts:", global_min_dist, global_max_dist
		print "Min and Max distance ranges:", min_range, max_range
		print "%i contacts above threshold" % (counter)
	
	return finalContacts

def getCutoffContacts(atoms, atomtypes, cutoff, chainDist, verbose=False, res2SS=None, selection=[], selectExclusive=False):
	contactAtoms = []
	
	for a in sorted(atoms.keys()):
		if atoms[ a ][ "aname" ] in atomtypes:
			contactAtoms.append(a)
	
	contacts = []

	for i in xrange(len(contactAtoms)):
		ai = contactAtoms[i]
		ri = atoms[ ai ][ "rnum" ]
		if selectExclusive and not ri in selection: continue
		if res2SS != None and res2SS[ri][0] == 'X':
			# print "skip ri",ri,res2SS[ri]
			continue

		for j in xrange(len(contactAtoms)):
			aj = contactAtoms[j]
			rj = atoms[ aj ][ "rnum" ]
			if selection != [] and selectExclusive and not rj in selection: continue
			
			if selection != [] and not selectExclusive and not ri in selection and not rj in selection: continue
			
			if res2SS != None and res2SS[rj][0] == 'X':
				# print "skip rj",rj,res2SS[rj]
				continue
			
			if ai < aj:
				if ri < rj - chainDist:
					d = atomDist(ai, aj, atoms)
					# print d
					if d <= cutoff:
						contacts.append([ai, aj, d])
						if verbose: print [ai, aj, d], ri, rj
	if verbose: print "Found %i CUTOFF contacts between %s atoms at chain separation=%i and cutoff distance = %f" % (len(contacts), str(atomtypes), chainDist, cutoff)
	
	return contacts

def getResidueDict(atoms):
	rd = {}
	
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in rd:
			rd[rnum].append(atoms[a])
		else:
			rd[rnum] = [atoms[a]]
	
	
	return rd

def getResidueContacts(atoms, cutoff, chainDist, verbose=False, res2SS=None, selection=[], selectExclusive=False):
	cutoff = cutoff
	chainsep = chainDist
	contacts = []
	residues = getResidueDict(atoms)
	n = len(residues)
	if verbose: print n, "residues"
	
	for i in xrange(n):
		if selectExclusive and not i + 1 in selection: continue
		ri = residues[i + 1]
		for j in xrange(n):
			if selection != [] and selectExclusive and not j + 1 in selection: continue
			
			if selection != [] and not selectExclusive and not i + 1 in selection and not j + 1 in selection: continue
			
			rj = residues[j + 1]
		
			if i < j - chainsep:
				mindist = 10000
				#mdpair = None
				cadist = None
				for ai in ri:
					if not (ai["aname"][0] == 'H' or (ai["aname"][0] in ['1', '2', '3', '4'] and ai["aname"][1] == 'H') or ai["aname"] in ["CA", "C", "O", "N"]): 
						# print ai.get_id()
						continue
					for aj in rj:
						if not (aj["aname"][0] == 'H' or (aj["aname"][0] in ['1', '2', '3', '4'] and aj["aname"][1] == 'H') or aj["aname"] in ["CA", "C", "O", "N"]): 
							continue
						d = atomDist(ai["anum"], aj["anum"], atoms)
					
						if cadist == None:
							if ai["aname"] == 'CA' and aj["aname"] == 'CA':
								cadist = d
						if d < mindist:
							mindist = d
							#mdpair = (ai["anum"], aj["anum"])
				if mindist <= cutoff:
					a1 = getAtomIDforResidueContact(i + 1, 'CA', atoms)
					a2 = getAtomIDforResidueContact(j + 1, 'CA', atoms)
					
					contacts.append([a1, a2, cadist])
					# print [i,j,mindist,cadist], mdpair[0], mdpair[1],a1,a2
	return contacts

def getResidueContactsSimple(pdbfile, chainsep):
	atoms = atomsFromPDB(pdbfile)
	residues = getResidueDict(atoms)
	n = len(residues)
	pairs = []
	dists = []
	for i in xrange(n):
		ri = residues[i + 1]
		for j in xrange(n):
			rj = residues[j + 1]
			if i < j - chainsep:
				mindist = 10000
				mdpair = None
				cadist = None
				natoms = 0
				for ai in ri:
					if not (ai["aname"][0] == 'H' or (ai["aname"][0] in ['1', '2', '3', '4'] and ai["aname"][1] == 'H') or ai["aname"] in ["CA", "C", "O", "N"]): 
						for aj in rj:
							if not (aj["aname"][0] == 'H' or (aj["aname"][0] in ['1', '2', '3', '4'] and aj["aname"][1] == 'H') or aj["aname"] in ["CA", "C", "O", "N"]): 
								
								d = atomDist(ai["anum"], aj["anum"], atoms)
					
								if cadist == None:
									if ai["aname"] == 'CA' and aj["aname"] == 'CA':
										cadist = d
								if d < mindist:
									mindist = d
									mdpair = (ai["anum"], aj["anum"])
									natoms += 1
							else:
								pass
					else:
						pass
				if  mindist == 10000:
					continue
				pairs.append((i + 1, j + 1))
				
				dists.append(mindist)
	return dists, pairs

def getAtomIDforResidueContact(resid, atype, atoms):
	for a in atoms:
		if atoms[a]["rnum"] == resid and atoms[a]["aname"] == atype:
			return a
	return None

def gro2pdbAtomMap(grofilename, pdbfilename, atomtypes):
	gro2pdb_heavy = {}
	groatoms = atomsFromGRO(grofilename)
	pdbatoms = atomsFromPDB(pdbfilename)
	assert (len(groatoms) - len(pdbatoms)) <= 1
	
	for i in groatoms.keys():
		resnr = groatoms[i]["resnr"]
		atype = groatoms[i]["atype"]
		if not atype in atomtypes: continue
		found = False
		for j in pdbatoms.keys():
			rnum = pdbatoms[j]["rnum"]
			aname = pdbatoms[j]["aname"]
			if atype == "OC1" and aname == "OXT":
				atype = "OXT"
			if resnr == rnum and atype == aname:
				gro2pdb_heavy[i] = j
				found = True
		if not found:
			print "could not find", i, resnr, atype
	return gro2pdb_heavy

def getShadowContacts(pdbfilename, outfilename, cutoff=6.0, chainDist=3, gmxpath="", scmpath="", res2SS=None, verbose=False, selection=[], selectExclusive=False):
	cleanupList = []
	
	grofilename = pdbfilename[:-4] + ".gro"
	topfilename = pdbfilename[:-4] + ".top"
	
	pdb2gmx = "%spdb2gmx -f %s -o %s -p %s -ff amber99sb-ildn -water tip3p -ignh -i tmpposre.itp" % (gmxpath, pdbfilename, grofilename, topfilename)
	if verbose: print pdb2gmx
	sub = subprocess.Popen(pdb2gmx, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	
	atoms = atomsFromGRO(grofilename)
	
	cleanupList.append(grofilename)
	cleanupList.append(topfilename)
	cleanupList.append("tmpposre.itp")
	
	scm = "java -jar %sSCM.jar -t %s -g %s -o %s --distance --coarse AACA -m shadow -s 1.0 -c %s --proteinDelta %s" % (scmpath, topfilename, grofilename, outfilename, str(cutoff), str(chainDist))
	if verbose: print scm
	sub = subprocess.Popen(scm, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	contacts = []
	
	gro2pdb = gro2pdbAtomMap(grofilename, pdbfilename, atomtypes=['CA'])
	
	for line in open(outfilename).readlines():
		tmp = line.strip().split()
		if len(tmp) == 5:
			a1 = int(tmp[1])
			a2 = int(tmp[3])
			d = float(tmp[4]) * 10.0
			r1 = atoms[ a1 ][ "resnr" ]
			r2 = atoms[ a2 ][ "resnr" ]

			if selection != [] and not selectExclusive and (not r1 in selection or not r2 in selection): continue
			
			if selection != [] and selectExclusive and not r1 in selection and not r2 in selection: continue
			
			a1 = gro2pdb[a1]
			a2 = gro2pdb[a2]
			if res2SS != None:
				if res2SS[r1][0] != 'X' and res2SS[r2][0] != 'X':
					contacts.append([ a1, a2, d ])
			else:
				contacts.append([ a1, a2, d ])
	
	cleanupList.append(outfilename)
	
	for i in cleanupList:
		if os.path.exists(i):
			os.remove(i)
	
	if verbose: print "Found %i SHADOW contacts in %s between CA atoms at chain separation=%i and cutoff distance = %f" % (len(contacts), pdbfilename, chainDist, cutoff)
	return contacts

def atomDist(a1, a2, atoms):
	d = math.sqrt(math.pow(atoms[a1]["x"] - atoms[a2]["x"], 2)
					+ math.pow(atoms[a1]["y"] - atoms[a2]["y"], 2)
					+ math.pow(atoms[a1]["z"] - atoms[a2]["z"], 2)
					)
	return d

def mixCutoffShadow(contacts_co, contacts_sh):
	contacts = []
	
	for c in contacts_co: contacts.append(c)

	for sh in xrange(len(contacts_sh)):
		a1_sh, a2_sh, d_sh = contacts_sh[sh]
		found = False
		for co in xrange(len(contacts_co)):
			a1_co, a2_co, d_co = contacts_co[co]
			
			if a1_sh == a1_co and a2_sh == a2_co:  # shadow contact is not in cutoff contacts
				found = True
				break
		if not found: 
			contacts.append(contacts_sh[sh])
	return contacts

def profasiAtomString(anum, atoms):
	rnum = atoms[anum]["rnum"]
	aname = atoms[anum]["aname"]	
	rname = atoms[anum]["rname"]
	return "0/%s/%s/_%s_" % (str(rnum - 1), rname, aname)

def profasiAtomStringNoResType(anum, atoms):
	rnum = atoms[anum]["rnum"]
	aname = atoms[anum]["aname"]	
	rname = atoms[anum]["rname"]
	return "0/%s//_%s_" % (str(rnum - 1), aname)

def atoms2atoms(atoms1, atoms2, atomtypes):
	a2a = {}
	for i in atoms1.keys():
		resnr = atoms1[i]["rnum"]
		atype = atoms1[i]["aname"]
		if not atype in atomtypes: continue
		found = False
		for j in atoms2.keys():
			rnum = atoms2[j]["rnum"]
			aname = atoms2[j]["aname"]
			
			if resnr == rnum and atype == aname:
				a2a[i] = j
				found = True
		if not found:
			print "could not find", i, resnr, atype
	return a2a

def getResNameFromNumber(atoms, rnum, aname):
	for a in atoms:
		if atoms[a]["rnum"] == rnum and atoms[a]["aname"]:
			return atoms[a]["rname"]
	return None

def getRes2AtomsDict(atoms):
	resdict = {}
	for a in atoms:
		rnum = atoms[a]["rnum"]
		if rnum in resdict:
			resdict[rnum].append(atoms[a])
		else:
			resdict[rnum] = [atoms[a]]
	return resdict

def writeProfasiContactsFMULTIGAUSS(contacts, atoms, atomtypes, outfilename, radius, steepness, width, depth, norm=0, label="", includeNonNative=False, chainSep=3, smoothGaussian=False, verbose=False):
	n = len(contacts)
	if norm > 0:
		depth = norm / float(n)
		if verbose: print "New well depth:", depth
	
	outfilename = outfilename + "_%icontacts_E%s" % (n, str(round(n * depth, 1))) + ".xml"
	
	out = open(outfilename, "w")
	
	if smoothGaussian:
		out.write("\n".join(['<sbm rid="%s">'%(label),
				'<formatted_data>'
				'  <format name=\"restraint\" type=\"$3\">'
				'  <atom1>$1</atom1>'
				'  <atom2>$2</atom2>'
				'  <parameters>'
				'    <low>$4</low>'
				'    <high>$5</high>'
				'    <radius>$6</radius>'
				'    <steepness>$7</steepness>'
				'    <width>$8</width>'
				'    <depth>$9</depth>'
				'  </parameters>'
				'  </format>'
				'  <data>'])+"\n")
	else:
		out.write("\n".join(['<sbm rid="%s">'%(label),
				'<formatted_data>',
				'  <format name=\"restraint\" type=\"$3\">',
				'  <atom1>$1</atom1>',
				'  <atom2>$2</atom2>',
				'  <parameters>',
				'    <minima>$4</minima>',
				'    <radius>$5</radius>',
				'    <steepness>$6</steepness>',
				'    <width>$7</width>',
				'    <depth>$8</depth>',
				'  </parameters>',
				'  </format>',
				'  <data>'])+"\n")
	
	seqdists = 0

	for ci in xrange(len(contacts)):
		c = contacts[ci]
		rnum1 = c[0]
		rnum2 = c[1]
		seqdists += abs(rnum1 - rnum2)
		atype1 = c[2]
		atype2 = c[3]
		dists = c[4]
		energies = [getEnergy(i, dists, radius, steepness, width, depth) for i in dists]
		assert all([i == -depth for i in energies])
		
		mindist = min(dists)
		if radius < 0:
			final_radius = mindist + radius
		else:
			final_radius = radius
		atomstring1 = "0/%s/%s/_%s_" % (str(rnum1 - 1), getResNameFromNumber(atoms, rnum1, atype1), atype1)  # profasiAtomString(a1, atoms)
		atomstring2 = "0/%s/%s/_%s_" % (str(rnum2 - 1), getResNameFromNumber(atoms, rnum2, atype2), atype2)  # profasiAtomString(a2, atoms)
		
		if smoothGaussian:
			out.write("    %s  %s  FMULTIGSMOOTH  %s  %s  %s %s %s %s\n" % (atomstring1, atomstring2, str(min(dists)), str(max(dists)), str(final_radius), str(steepness), str(width), str(depth)))
		else:	
			out.write("    %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n" % (atomstring1, atomstring2, ",".join([str(i) for i in dists]), str(final_radius), str(steepness), str(width), str(depth)))
	
	N = getLength(atoms)
	if verbose: print "CO = ", seqdists / float(n * N), seqdists, n, N
	
	if includeNonNative:
		nnCounter = 0
		nCounter = 0
		resdict = getRes2AtomsDict(atoms)
		
		for i in xrange(N):
			ai = None
			for a in resdict[i + 1]: 
				if a["aname"] in atomtypes: ai = a
			assert ai != None
			
			for j in xrange(N):
				assert isinstance(i , (int))
				if i < j - chainSep:
					aj = None
					for a in resdict[j + 1]: 
						if a["aname"] in atomtypes: aj = a
					assert aj != None
					
					isNative = False
					for ci in xrange(len(contacts)):
						c = contacts[ci]
						rnum1 = c[0]
						rnum2 = c[1]
						atype1 = c[2]
						atype2 = c[3]
						if rnum1 == i + 1 and rnum2 == j + 1 and ai["aname"] == atype1 and aj["aname"] == atype2:
							isNative = True
					
					if not isNative:
						nnCounter += 1
						atomstring1 = "0/%s/%s/_%s_" % (str(i), ai["rname"], ai["aname"])
						atomstring2 = "0/%s/%s/_%s_" % (str(j), aj["rname"], aj["aname"])
		
						out.write("  %s  %s  FMULTIGAUSS  %s  %s  %s %s %s\n" % (atomstring1, atomstring2, ",".join(["0" for x in dists]), "4", "1", "0", "0"))
					else:
						nCounter += 1
		if verbose: 
			print nnCounter, "non-native contacts written! (chain sep = %i)" % chainSep
			print nCounter, "native contacts"
		assert nCounter == n
		assert nnCounter + nCounter == (N * N) / 2.0 - N * 0.5 - sum([N - xx for xx in range(1, chainSep + 1)]), str(nnCounter + nCounter) + " " + str((N * N) / 2.0 - N * 0.5 - sum([N - xx for xx in range(1, chainSep + 1)]))
	out.write("\n".join([
					'  </data>',
					'</formatted_data>',
					'</sbm>']))
	if verbose: print "Expected energy:", n * depth
	
	return outfilename

def writeChimeraPseudobonds(contacts, outfilename, legend=None, verbose=False):

	outfile = open(outfilename, 'w')
	
	if legend != None: assert len(contacts) == len(legend), str(len(contacts)) + " " + str(len(legend))
	
	maxrange = 0
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		dists = c[4]
		rn = max(dists) - min(dists)
		if rn > maxrange:
			maxrange = rn
	#print "maxrange", maxrange
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		rnum1 = c[0]
		rnum2 = c[1]
		atype1 = c[2]
		atype2 = c[3]
		dists = c[4]
		rn = float(max(dists) - min(dists))
		
		label = str(round(min(dists), 2)) + "-" + str(round(max(dists), 2))
		
		if legend == None:
			# color = "blue"
			if maxrange > 0:
				gray = int(math.ceil((1.0 - (rn / maxrange)) * 65535))
			else:
				gray = 65535
			# print gray
			hx = hex(gray)
			color = "#" + str(hx).replace("0x", "") * 3
			# color = color.upper()
			# print gray,color,hx
		else:
			lstring = legend[ci]
			if lstring == "s1":
				color = "blue"
			elif lstring == "s2":
				color = "red"
			else:
				color = "magenta"
		pseudobond = "#0:%s@%s #0:%s@%s %s %s" % (str(rnum1), atype1.lower(), str(rnum2), atype2.lower(), color, label)
		outfile.write(pseudobond + "\n")

	outfile.close()

def getEnergy(x, minima, radius, steepness, width, depth):
	
	wells = 1.0
	for m in minima:
		if width != 0.0:
			Gij = -math.exp(-math.pow(x - m, 2.0) / (2 * width * width))
	
		else:
			Gij = 0.0
	
		wells *= (1 + Gij)

	Rij = steepness * math.pow(radius / x, 12.0)
	repulsion = (1 + (Rij / depth))
		
	return (depth * repulsion * wells) - depth

def plotContactMap(n, contacts, inPDBFilename, cutoff, chainSep, atoms, contactType, fname="contactmap.png", verbose=False):
	cdict = {}
	for c in contacts:
		r1 = c[0]
		r2 = c[1]
		dists = c[4]
		
		if len(dists) == 1:
			cdict[(r1, r2)] = 1.0
		else:
			cdict[(r1, r2)] = max(dists) - min(dists)  # plot RANGE of CA distances
	
	matrix = [["NaN" for j in xrange(n + 1)] for i in xrange(n + 1) ]
	
	for i in xrange(n + 1):
		for j in xrange(n + 1):
			if (i, j) in cdict:
				matrix[i][j] = cdict[(i, j)]
				matrix[j][i] = cdict[(i, j)]
	
	open(fname + ".data", "w").write("\n".join([",".join([str(j) for j in i[1:]]) for i in matrix[1:]]))
	
	matrix = [[matrix[i][j] if matrix[i][j] != "NaN" else 0.0 for j in xrange(n + 1)] for i in xrange(n + 1) ]
	
	title = "%s; cutoff=%s; csep=%i; n=%i\n%s" % (inPDBFilename, str(round(cutoff, 2)), chainSep, len(contacts), contactType)

	imax = p.matshow(matrix, cmap=cm.afmhot, aspect='equal', origin='lower')
	
	p.plot([1, n], [1, n], color='k')
	p.title(title)
	cbar = p.colorbar()
	p.xlabel("residue position")
	p.ylabel("residue position")
	p.xlim(0.5, n + 0.5)
	p.ylim(0.5, n + 0.5)
	p.grid(color="0.5", linestyle=':', linewidth=1.0)
	cbar.set_label("CA dist. range")
	p.savefig(fname, dpi=600)
	# p.show()
	if verbose: print "wrote contact map PNG:", fname

def getLength(atoms):
	maxrnum = 0
	for a in atoms:
		if atoms[a]["rnum"] > maxrnum:
			maxrnum = atoms[a]["rnum"]
	return maxrnum

def getPrfSCcontactEnergy(tconf, pdb, selection, parameters, prf_loc="~/workspace/PROFASI/app/bin/"):
	# ~/workspace/PROFASI/app/bin/prf_energies 1PGA_min_rmsd.tconf --box_length 1000 --skip_energy FF08:ExVol,FF08:Bias,FF08:LocExVol,FF08:HBMM,FF08:TorsionTerm --add_chain_pdb 1 1PGA_min_rmsd.pdb::A,45,45:A,49,49 -force_field FF14 --cationpi 10
	cmd = "%sprf_energies %s --box_length 1000 --skip_energy FF08:ExVol,FF08:Bias,FF08:LocExVol,FF08:HBMM,FF08:TorsionTerm --add_chain_pdb 1 %s%s -force_field FF14 %s" % (prf_loc, tconf, pdb, selection, parameters)
	sub = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	sub.wait()
	out = sub.stdout.read()
	err = sub.stderr.read()
	print err
	print out
	HBMS = None
	Hydrophobicity = None
	ChargedSCInteraction = None
	CationPi = None
	Total = None
	start = False
	for line in out.split("\n"):
		if ">>>>>>>>>>" in line:
			# print "start"
			start = True
			continue
		
		if start:
			tmp = line.strip().split()
			if len(tmp) == 3:
				# print tmp
				if tmp[0] == "HBMS": HBMS = float(tmp[2])
				elif tmp[0] == "Hydrophobicity": Hydrophobicity = float(tmp[2])
				elif tmp[0] == "ChargedSCInteraction": ChargedSCInteraction = float(tmp[2])
				elif tmp[0] == "CationPi": CationPi = float(tmp[2])
				elif tmp[0] == "Total": Total = float(tmp[2])
	return (HBMS, Hydrophobicity, ChargedSCInteraction, CationPi, Total)

def getSCenergies4Contacts(contacts, pdb, atoms):
	energies = []
	print "Profasi SC energies per contact:"
	for c in contacts:
		
		r1 = c[0]
		r2 = c[1]
		a1 = c[2]
		a2 = c[3]
		selection = "::A,%i,%i:A,%i,%i" % (r1, r1, r2, r2)
		parameters = "--cationpi 10"
		tconf = pdb[:-4] + ".tconf"
		sc_energy = getPrfSCcontactEnergy(tconf, pdb, selection, parameters, prf_loc="~/workspace/PROFASI/app/bin/")
		print c, sc_energy, getResNameFromNumber(atoms, r1, a1), getResNameFromNumber(atoms, r2, a2)
		energies.append(sc_energy)
	return energies

def multiGaussSmoothFunc(x, dists, eps, w):
	vals = []
	low = min(dists)
	high = max(dists)
	for i in x:
		if i <= low:
			vals.append(eps * (1.0 - numpy.exp(-(numpy.power(i - low, 2.0) / (2.0 * w * w)))) - eps)
		elif i >= high:
			vals.append(eps * (1.0 - numpy.exp(-(numpy.power(i - high, 2.0) / (2.0 * w * w)))) - eps)
		else:
			vals.append(-eps)
	return numpy.array(vals)

def multiGaussFunc(x, dists, eps, w):
	prod = eps
	for i in dists:
		prod *= 1.0 - numpy.exp(-(numpy.power(x - i, 2.0) / (2.0 * w * w))) 
	return eps * prod - eps

def LennardJones(x, dist, eps):
	return eps * (numpy.power(dist / x, 12) - 2 * numpy.power(dist / x, 6))

def multiLJ(x, dists, eps):
	return sum([LennardJones(x, d, eps) for d in dists])

def plotMultiGaussianContacts(contacts, eps, w, fname, smoothGaussian=False, verbose=False):
	l = len(contacts)
	fig = p.Figure(figsize=(25, 25), dpi=800)
	p.suptitle("")
	sq = int(math.ceil(math.sqrt(len(contacts))))
	
	globmin = 10e10
	globmax = -10e10
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		r1 = c[0]
		r2 = c[1]
		dists = c[4]
		if min(dists) < globmin: globmin = min(dists)
		if max(dists) > globmax: globmax = max(dists)
	
	for ci in xrange(len(contacts)):
		c = contacts[ci]
		r1 = c[0]
		r2 = c[1]
		dists = c[4]
		
		ax = p.subplot(sq, sq, ci + 1)
		
		x = numpy.arange(math.floor(globmin) - 1.0, math.ceil(globmax) + 1.0, 0.01)
		ax.set_ylim(-eps - 0.1, 0.2)
		ax.set_xlim(math.floor(globmin) - 1.0, math.ceil(globmax) + 1.0)
		
		textsize = 5
		
		if ci in [i for i in xrange(len(contacts)) if i == 0 or i % sq == 0 ]: 
			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(textsize - 1)
			ax.set_ylabel("epsilon", size=textsize)
		else: 
			ax.set_yticklabels([])
		
		if ci in [i for i in range(len(contacts) - sq, len(contacts)) ]:	
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(textsize - 1)
			ax.set_xlabel("dij", size=textsize)
		else: 
			ax.set_xticklabels([])
		
		for d in dists:
			p.plot(x, LennardJones(x, d, eps), linestyle='-', linewidth=0.1, c='red')
		p.plot(x, multiLJ(x, dists, eps) / len(dists), linestyle='-', linewidth=0.4, c='green')
		
		if smoothGaussian:
			p.plot(x, multiGaussSmoothFunc(x, dists, eps, w), linewidth=0.8, c='blue')
		else:
			p.plot(x, multiGaussFunc(x, dists, eps, w), linewidth=0.8, c='blue')
		
		p.text(globmin, 0.3, "%i: (%i,%i)" % (ci + 1, r1, r2), size=textsize)
	
	p.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, wspace=0.1, hspace=0.4)
	p.savefig(fname, dpi=600)
	# p.show()
	
def createLinkSubfoldersKmeans(pdbfolder, centroids, feature, clusterMembers, shortlinks=False):
	assert len(clusterMembers) == len(centroids)
	
	rmExistingClustSubFolders(pdbfolder, "k%imeans_%s" % (len(centroids), feature))
	
	options = ""
	centroidFile = open(pdbfolder + "".join(options.replace("-", "_").split()) + "_k%imeans_%s_centroids.txt" % (len(centroids), feature), 'w')
	centroids_subfolder = os.path.join(pdbfolder, "k%imeanscentroids_%s" % (len(centroids), feature))
	
	if os.path.exists(centroids_subfolder): 
		shutil.rmtree(centroids_subfolder)
	
	os.mkdir(centroids_subfolder)
	
	subfolders = []
	
	for c in xrange(len(clusterMembers)):
		
		centroidFile.write(centroids[c] + "\n")
		
		if shortlinks:
			cenlink = os.path.abspath(os.path.join(centroids_subfolder, "cen%s.pdb" % str(c).zfill(len(str(len(centroids))))))
		else:
			cenlink = os.path.abspath(os.path.join(centroids_subfolder, os.path.basename(centroids[c])))
		
		censrc = os.path.abspath(centroids[c])
		
		if not os.path.exists(cenlink):
			if os.path.lexists(cenlink):
				os.remove(cenlink)
			os.symlink(censrc, cenlink)
		
		cluster_subfolder = os.path.join(pdbfolder, "k%imeans_%s%s" % (len(centroids), feature, str(c).zfill(len(str(len(clusterMembers))))))
		if not os.path.exists(cluster_subfolder): 
			os.mkdir(cluster_subfolder)
		subfolders.append(cluster_subfolder)
		for mi in xrange(len(clusterMembers[c])):
			m = clusterMembers[c][mi]
			if shortlinks:
				memlink = os.path.abspath(os.path.join(cluster_subfolder, "mem%s.pdb" % str(mi).zfill(len(str(len(clusterMembers[c]))))))
			else:
				memlink = os.path.abspath(os.path.join(cluster_subfolder, os.path.basename(m)))
			memsrc = os.path.abspath(m)
			# link_legend_file.write("%s %s\n"%(memlink, memsrc))
			if not os.path.exists(memlink):
				if os.path.lexists(memlink):
					os.remove(memlink)
				os.symlink(memsrc, memlink)
	return subfolders, centroids_subfolder


def createLinkSubfoldersMaxclust(pdbfolder, options, centroids, clusterMembers, id2pdb, shortlinks=False):
	centroidFile = open(pdbfolder+"".join(options.replace("-","_").split())+"_centroids.txt",'w')
	centroids_subfolder = os.path.join(pdbfolder,"centroids")


	if os.path.exists(centroids_subfolder): 
		shutil.rmtree(centroids_subfolder)

	os.mkdir(centroids_subfolder)

	rmExistingClustSubFolders(pdbfolder)

	link_legend_file = open(pdbfolder+"".join(options.replace("-","_").split())+"_linklegend.txt",'w')
	
	# write "cluster" 0 - contains unassigned pdbs
	
	
	for c in xrange(len(centroids)+1):
		if c >0:
			centroidFile.write(centroids[c][3]+"\n")
			
			if shortlinks:
				cenlink = os.path.abspath(os.path.join(centroids_subfolder, "cen%s.pdb"%str(c).zfill(len(str(len(centroids))))))
			else:
				cenlink = os.path.abspath(os.path.join(centroids_subfolder, os.path.basename(centroids[c][3])))
			censrc = os.path.abspath(centroids[c][3])
			link_legend_file.write("%s %s\n"%(cenlink, censrc))
			if not os.path.exists(cenlink):
				if os.path.lexists(cenlink):
					os.remove(cenlink)
				os.symlink(censrc, cenlink)
	
		cluster_subfolder = os.path.join(pdbfolder,"clust%s"%str(c).zfill(len(str(len(centroids)))))
		if not os.path.exists(cluster_subfolder): 
			os.mkdir(cluster_subfolder)
	
		for m in clusterMembers[c]:
			if shortlink:
				memlink = os.path.abspath(os.path.join(cluster_subfolder, "mem%s.pdb"%str(m).zfill(len(str(len(id2pdb))))))
			else:
				memlink = os.path.abspath(os.path.join(cluster_subfolder, os.path.basname(id2pdb[m])))
			memsrc = os.path.abspath(id2pdb[m])
			link_legend_file.write("%s %s\n"%(memlink, memsrc))
			if not os.path.exists(memlink):
				if os.path.lexists(memlink):
					os.remove(memlink)
				os.symlink(memsrc, memlink)
def filterSymLinksBySource(flist, filtstr):
	filenames = []
	
	for f in flist:
		if os.path.islink(f):
			if filtstr in os.readlink(f):
				filenames.append(f)
	
	return filenames

def addFile2Data(pdb2data, datafilename):
	df = open(datafilename)
	h = df.readline().strip().split()
	for line in df.readlines():
		tmp = line.strip().split(",")
		assert not tmp[0] in pdb2data
		pdb2data[tmp[0]] = [ float(i) for i in tmp[1:] ]
	return h

def rmExistingClustSubFolders(folder, prefix):
	flist = glob.glob(os.path.join(folder, "%s*" % prefix))
	for f in flist:
		if os.path.isdir(f):
			shutil.rmtree(f)	

def doKmeans4Feature(infolder, k, feature, njobs=2, alwaysOverwrite=False, verbose=False, ref_pdb=None, atom_types=None):
	from msmbuilder.cluster import KMeans
	from msmbuilder.dataset import dataset
	from msmbuilder.featurizer import ContactFeaturizer, DihedralFeaturizer, SuperposeFeaturizer, RMSDFeaturizer
	import mdtraj as md
	
	assert feature in ["superdihedralcontacts", "contacts", "dihedrals", "superposition", "rmsd"]
	
	resultsfilename = "%s_kmeans%i_%s.txt" % (infolder.replace("/", "__"), k, feature)
	
	if not os.path.exists(resultsfilename):
	
		filelist = glob.glob(os.path.join(infolder, "*.pdb"))

		ds = []
	
		#cf = ContactFeaturizer()
		df = DihedralFeaturizer()
		sf = None
		rf = None
		r = None
		if feature in ["superposition", "rmsd", "superdihedralcontacts"]: 
			assert ref_pdb != None
			
			if atom_types == None:
				iatoms = [None for i in ref_pdb]
			else:
				iatoms = [getAtomList(i, atom_types) for i in ref_pdb]
			
			if feature in ["superposition", "superdihedralcontacts"]:
				sf = []
				for i  in xrange(len(ref_pdb)):
					r = md.load(ref_pdb[i])
				sf.append(SuperposeFeaturizer(iatoms[i], r))
			else:
				r = md.load(ref_pdb)
				rf = RMSDFeaturizer(r)
		counter = 0
		for f in filelist:
			sys.stdout.write("\r%s%%" % str(round(float(100.0 * counter) / len(filelist), 2)))
			sys.stdout.flush()
			t = md.load(f)
			
			if feature == "contacts":
				feats, pairs = getResidueContacts(f, 3)
				feats = [feats]
			elif feature == "dihedrals":
				feats = df.partial_transform(t)
			elif feature == "superposition":
				feats = []
				for s in sf:
					feats.extend(s.partial_transform(t)[0])
			elif feature == "rmsd":
				feats = rf.partial_transform(t)
			elif feature == "superdihedralcontacts":
				feats = []
				for s in sf:
					feats.extend(s.partial_transform(t)[0])
				
				feats.extend(df.partial_transform(t)[0])
				confeats, pairs = getResidueContacts(f, 3)
				feats.extend(confeats)
			
			ds.append(numpy.array([feats]))
			counter += 1
	
		if verbose: print "dataset dims:", len(ds), len(ds[0]), len(ds[0][0])

		cluster = KMeans(n_clusters=k, n_jobs=njobs)
		cluster.fit(ds)
	
		centers = cluster.cluster_centers_
		inertia = cluster.inertia_
	
		if verbose: 
			print "centers:", len(centers), len(centers[0])
			print "inertia:", inertia
	
		# process clusters

		pdbclusters = [[] for i in xrange(k)]

		for l in xrange(len(cluster.labels_)):
			pdbclusters[cluster.labels_[l][0]].append(filelist[l])
	
		centroids = []

		for i in xrange(len(centers)):
			print "center", i
			nearest = None
			nearest_dist = 10e10
			for j in xrange(len(ds)):
				x = ds[j]
				dist = numpy.linalg.norm(x - centers[i]) 
				d2 = dist * dist
				if d2 < nearest_dist:
					nearest_dist = d2
					nearest = j
			print nearest, filelist[nearest]
			centroids.append(filelist[nearest])

		print centroids
		
		r = open(resultsfilename, 'w')
		r.write(",".join(centroids) + "\n")
		for i in xrange(len(pdbclusters)):
			r.write(",".join(pdbclusters[i]) + "\n")
		r.close()
	else:
		pdbclusters = []
		r = open(resultsfilename)
		print "Loading", resultsfilename
		centroids = r.readline().strip().split(",")
		for line in r.readlines():
			tmp = line.strip().split(",")
			pdbclusters.append(tmp)
		assert len(pdbclusters) == k
	return pdbclusters, centroids

def getAdjacencyMatrixFromFeature(pdblist, cutoff, feature, ref_pdb=[], atom_types=None):
	from msmbuilder.cluster import KMeans
	from msmbuilder.dataset import dataset
	from msmbuilder.featurizer import ContactFeaturizer, DihedralFeaturizer, SuperposeFeaturizer, RMSDFeaturizer
	import mdtraj as md
	
	dlist = []
	
	for p in pdblist:
		if feature == "contacts":
			feat, pairs = getResidueContacts(p, 3)
		elif feature == "dihedrals":
			df = DihedralFeaturizer()
			t = md.load(p)
			feat = df.partial_transform(t)[0]
		elif feature == "superposition":
			if atom_types == None:
				iatoms = [None for i in ref_pdb]
			else:
				iatoms = [getAtomList(i, atom_types) for i in ref_pdb]
			
			feat = []
			for i  in xrange(len(ref_pdb)):
				r = md.load(ref_pdb[i])
				sf = SuperposeFeaturizer(iatoms[i], r)
				t = md.load(p)
				feat.extend(sf.partial_transform(t)[0])
		dlist.append(feat)
	
	A = [[None for j in xrange(len(pdblist))] for i in xrange(len(pdblist)) ]
	print len(pdblist), len(dlist), len(dlist[0])
	allrms = []
	for i in xrange(len(pdblist)):
		
		for j in xrange(len(pdblist)):
			if i != j:
				squares = []
				for d in xrange(len(dlist[i])):
					squares.append(math.pow(dlist[i][d] - dlist[j][d], 2))
				rms = math.sqrt(sum(squares) / float(len(dlist[0])))
				
				if rms <= cutoff:
					allrms.append(rms)
					A[i][j] = rms
	
	maxrms = float(max(allrms))
	minrms = float(min(allrms))
	print "mean/min/max rms, n:", numpy.mean(allrms), minrms, maxrms, len(allrms)
	
	allrms2 = []
	for i in xrange(len(pdblist)):
		
		for j in xrange(len(pdblist)):
			if A[i][j] != None:
				rms2 = 1.0 - ((A[i][j] - minrms) / (maxrms - minrms))
				# print A[i][j], rms2
				A[i][j] = rms2
				allrms2.append(rms2)
			else:
				A[i][j] = 0.0
			
	
	print "mean/min/max rms, n:", numpy.mean(allrms2), min(allrms2), max(allrms2)
	
	return A

def splitModels(pdbfilename, verbose=False):
	hasModels = False
	models = []
	model_nr = 0
	atomlines = []
	for line in open(pdbfilename).readlines():
		linetype = line[:6]
		
		if linetype == "MODEL ":
			hasModels = True
			model_nr = int(line.replace("MODEL", ""))
			if verbose: print "model nr ", model_nr
			atomlines = []
		elif linetype == "ENDMDL":
			assert model_nr != 0, model_nr
			modelfilename = pdbfilename[:-4] + "_m%i.pdb" % model_nr
			models.append(modelfilename)
			
			modelfile = open(modelfilename, "w")

			assert len(atomlines) != 0
			for l in atomlines:
				modelfile.write(l)
			modelfile.write("END   ")
			modelfile.close()
		elif linetype == "ATOM  ":
			atomlines.append(line)
		else:
			pass
	if hasModels:
		return models  # returns a list of the new filenames, one model per file
	else:
		return [pdbfilename]

def joinModels(listOfFilenames, outfilename, cleanUp=False, verbose=False):
	if verbose: print listOfFilenames
	outfile = open(outfilename, 'w')
	counter = 1
	for infilename in listOfFilenames:
		if not os.path.exists(infilename):
			print "WARNING: file does not exist! ", infilename
			continue
		infile = open(infilename)

		spaces = (9 - int(math.log10(counter))) * " "
		outfile.write("MODEL" + spaces + str(counter) + "\n")
		outfile.write("REMARK original_filename %s\n" % (infilename))
		for line in infile.readlines():

			if line[:4] == "ATOM":
				outfile.write(line)

		outfile.write("ENDMDL                                                                          \n")
		counter += 1
		infile.close()
	outfile.write("END                                                                             \n")
	outfile.close()
	if cleanUp: # CAUTION: removes all PDB files in list!
		for l in listOfFilenames:
			if os.path.exists(l): os.remove(l)

def getAtomList(pdb, atom_types):
	filtered = []
	atoms = atomsFromPDB(pdb)
	for a in atoms:
		if atoms[a]["aname"].strip() in atom_types:
			filtered.append(int(atoms[a]["anum"]))
	return filtered

def parseATOMline(line, returnDict=True): 
	d = {}
	d["ltype"] = line[:6]
	try:
		d["anum"] = int(line[6:11])
		d["aname"] = line[12:16].strip()
		d["aname_exact"] = line[12:16]
		d["altloc"] = line[16]
		d["rname"] = line[17:20].strip()
		d["chain"] = line[21]
		d["rnum"] = int(line[22:26])
		d["insert"] = line[26]
		if d["ltype"] in ["ATOM  ", "HETATM"]:
			d["x"] = float(line[30:38])
			d["y"] = float(line[38:46])
			d["z"] = float(line[46:54])
			d["occupancy"] = float(line[54:60])
			d["bfactor"] = float(line[60:66])
			d["segment"] = line[72:76].strip()
			d["element"] = line[76:78].strip()
			d["charge"] = line[78:80].strip()
	except ValueError as ve:
		print ve
		print line
		sys.exit(1)
	if returnDict:
		return d
	else:
		return (d["ltype"], d["anum"], d["aname"], d["altloc"], d["rname"], d["chain"], d["rnum"], d["insert"], d["x"], d["y"], d["z"], d["occupancy"], d["bfactor"], d["segment"], d["element"], d["charge"])

def writeATOMdict2Line(d):
	ltype = d["ltype"] 
	line = ltype
	line += "{:>5}".format(str(d["anum"]))
	line += " "
	line += "{:<4}".format(d["aname_exact"])
	line += "{}".format(d["altloc"])
	line += "{}".format(d["rname"])
	line += " "
	line += "{}".format(d["chain"])
	line += "{:>4}".format(d["rnum"])
	line += "{}".format(d["insert"])
	if ltype in ["ATOM  ", "HETATM"]:
		line += "   "
		line += "{:>8.3f}".format(d["x"])
		line += "{:>8.3f}".format(d["y"])
		line += "{:>8.3f}".format(d["z"])
		line += "{:>6.2f}".format(d["occupancy"])
		line += "{:>6.2f}".format(d["bfactor"])
		line += "     "
		line += "{:<4}".format(d["segment"])
		line += " "
		line += "{:>2}".format(d["element"])
		line += "{:>2}".format(d["charge"])
	elif ltype == "TER   ":
		line += "                                                      "
	return line + "\n"

def prf_convert(infile, outfile, prf_loc):
	exec_string = "%sprf_convert %s -o %s" % (prf_loc, infile, outfile)
	sub = subprocess.Popen(exec_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out = sub.stdout.read()
	err = sub.stderr.read()
	if err.strip() != "": 
		print out
		print err

def regularize(inpdb, prf_loc, options="", verbose=False):
	
	exec_string = "%sregularize %s %s" % (prf_loc, inpdb, options)
	sub = subprocess.Popen(exec_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out = sub.stdout.read()
	err = sub.stderr.read()
	
	if verbose:#if err.strip() != "": 
		print out
		print err
	
	# PROFASI structure with smallest RMSD from original
	prf_convert("min_rmsd.xml", inpdb[:-4] + "_min_rmsd.pdb", prf_loc)
	prf_convert("min_rmsd.xml", inpdb[:-4] + "_min_rmsd.tconf", prf_loc)
	shutil.move("min_rmsd.xml", inpdb[:-4] + "_min_rmsd.xml")
	
	# PROFASI structure with lowest energy'
	prf_convert("min_etot.xml", inpdb[:-4] + "_min_etot.pdb", prf_loc)
	prf_convert("min_etot.xml", inpdb[:-4] + "_min_etot.tconf", prf_loc)
	shutil.move("min_etot.xml", inpdb[:-4] + "_min_etot.xml")

def splitModelsAndChains(pdbfilename, chainlabel="A"):
	hasModels = False
	models = []
	model_nr = 0
	chains = []
	currchain = None
	atomlines = []
	atomcounter = 0
	for line in open(pdbfilename).readlines():
		
		linetype = line[:6]
		
		if linetype == "MODEL ":
			hasModels = True
			model_nr = int(line.replace("MODEL", ""))
			print "model nr ", model_nr
			atomlines = []
			atomcounter = 0
			chains = []
			currchain = None
		elif linetype == "ENDMDL":
			assert model_nr != 0, model_nr
			
			print len(atomlines), "atoms found in chain", currchain
			chains.append(atomlines)
			
			print len(chains), "chains found"
			for ci in xrange(len(chains)):
				modelfilename = pdbfilename[:-4] + "_m%i_c%i.pdb" % (model_nr, ci + 1)
				models.append(modelfilename)
				modelfile = open(modelfilename, "w")

				assert len(chains[ci]) != 0
				for l in chains[ci]:
					modelfile.write(l)
				modelfile.write("END   ")
				modelfile.close()
			chains = []
		elif linetype in ["ATOM  ", "TER   ", "HETATM"]:
			a = parseATOMline(line)
			c = a["chain"]
			# print c
			if currchain == None:
				currchain = c
			else:
				if c != currchain:
					print "new chain:", c, " old chain:", currchain
					chains.append(atomlines)
					print len(atomlines), "atoms found in chain", currchain
					atomlines = []
					atomcounter = 0
					currchain = c
			a["chain"] = chainlabel
			a["anum"] = atomcounter + 1
			atomlines.append(writeATOMdict2Line(a))
			atomcounter += 1
		else:
			pass
	
	print len(models), "models found in", pdbfilename
	
	if not hasModels and len(models) == 0:
		print len(chains), "chains found"
		print len(atomlines), "atoms found in chain", currchain
		chains.append(atomlines)
		
		print len(chains), "chains found"
		for ci in xrange(len(chains)):
			modelfilename = pdbfilename[:-4] + "_m%i_c%i.pdb" % (1, ci + 1)
			models.append(modelfilename)
			modelfile = open(modelfilename, "w")

			assert len(chains[ci]) != 0
			for l in chains[ci]:
				modelfile.write(l)
			modelfile.write("END   ")
			modelfile.close()
	
	return models  # returns a list of the new filenames, one chain and model per file

def downloadPDB(code):
	import urllib
	filename = "%s.pdb"%(code.upper())
	urllib.urlretrieve ("http://www.rcsb.org/pdb/files/%s"%filename, filename)
	
	if os.path.exists(filename):
		return True
	else:
		return False

def getList4Folder(folder, options="", n=None):
	assert os.path.isdir(folder)
	listfile = folder+"".join(options.replace("-","_").split())+"_maxcluster_list.txt"
	flist = glob.glob(folder+"/*.pdb")
	if n!=None and len(flist)>n:
		flist = random.sample(flist, n)
		
	lf = open(listfile,"w")
	lf.write("\n".join(flist)+"\n")
	lf.close()
	
	return listfile, len(flist)

# run maxcluster on specified folder and save output in text file
def maxcluster(folder, listfile, options="", mcexepath="maxcluster"):
	exe = "%s -l %s %s"%(mcexepath, listfile, options)
	#else:
	#	exe = "cd %s; %s -l %s %s; cd .."%(folder, mcexepath, listfile, options)
	print exe
	p = subprocess.Popen(exe, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out = p.stdout.read()
	err = p.stdout.read()
	
	if err.strip() != "":
		print "Maxcluster error:"
		print err
		return None
	
	outfilename = folder + "".join(options.replace("-","_").split()) + "_maxcluster_out.txt"
	o = open(outfilename,'w')
	o.write(out)
	o.close()
	
	return outfilename

# process text output
#returns:
#centroids --> clusterNo:(centroidNo, clusterSize, clusterSpread, centroidPdb)
#clusterMembers --> clusterNo:[pdbId]
#id2pdb --> pdbId:pdb
def parseMaxclusterOutput(filename, options="",labels=[]):
	lines = open(filename).readlines()
	centroids_start_index = None
	centroids_end_index = None
	for l in xrange(len(lines)):
		if lines[l].strip()=="INFO  : Centroids" and lines[l+1].strip()=="INFO  : ======================================":
			assert lines[l-2].strip()=="INFO  : Nearest Neighbour clustering", lines[l-2] # TODO: extend for output from other clustering methods
			assert lines[l+2].strip()=="INFO  : Cluster  Centroid  Size        Spread", lines[l+2]
			centroids_start_index = l+3
			break
	
	#find cluster properties and centroids
	centroids={} # clusterNo:(centroidNo, clusterSize, clusterSpread, centroidPdb)
	cluster_info = None
	for l in range(centroids_start_index, len(lines)):
		if centroids_end_index == None and lines[l].strip() != "INFO  : ======================================":
			tmp = lines[l].strip().split()
			clusterNo = int(tmp[2])
			centroidNo = int(tmp[4])
			clusterSize = int(tmp[5])
			clusterSpread = float(tmp[6])
			centroidPdb = tmp[7]
			centroids[clusterNo] = (centroidNo, clusterSize, clusterSpread, centroidPdb)
			print clusterNo, centroids[clusterNo]
		else:
			centroids_end_index = l-1
			cluster_info = lines[l+1].strip()
			cis = cluster_info.split()
			nclust = int(cis[2])
			assert nclust==len(centroids), cluster_info
			assert cis[5]=="Threshold", cis[5]
			threshold = float(cis[6])
			assigned = int(cis[8])
			total = int(cis[10].strip(")"))
			break
	print nclust, threshold, assigned, total
	
	histo = {}
	
	clusterMembers = {}
	id2pdb = {}
	for l in range(centroids_end_index+5, len(lines)):
		tmp = lines[l].strip().split()
		if len(tmp) == 6:
			pdbNo = int(tmp[2])
			clusterNo = int(tmp[4])
			pdb = tmp[5]
			assert not pdbNo in id2pdb
			id2pdb[pdbNo] = pdb
			if labels!=[]:
				assert sum([1 for i in xrange(len(labels)) if labels[i] in pdb])==1, pdb
				if clusterNo in histo:
					
					for i in xrange(len(labels)):
						if labels[i] in pdb:
							histo[clusterNo][i] += 1
				else:
					histo[clusterNo] = [0 for i in labels]
					for i in xrange(len(labels)):
						if labels[i] in pdb:
							histo[clusterNo][i] += 1
			if clusterNo in clusterMembers:
				clusterMembers[clusterNo].append(pdbNo)
			else:
				clusterMembers[clusterNo] = [pdbNo]
	
	assert len(id2pdb)==total
	for c in clusterMembers:
		if c!=0:
			assert len(clusterMembers[c])==centroids[c][1]
			assert centroids[c][0] in id2pdb
		else:
			assert len(clusterMembers[c])==total-assigned
		#print c, len(clusterMembers[c])
	return centroids, clusterMembers, id2pdb, histo

def rmExistingClustSubFolders(folder):
	flist = glob.glob(os.path.join(folder, "clust*"))
	for f in flist:
		if os.path.isdir(f):
			shutil.rmtree(f)


def mergeStrings(slist):
	s=""
	for i in xrange(min( [len(sj) for sj in slist] )):
		if all( [slist[0][i]==slist[j][i] for j in xrange(len(slist))] ):
			s += slist[0][i]
		else:
			s += "x"
	return s
def mergeFoldersLinks(folders,n=None):
	newname = mergeStrings(folders)+"_merged"
	assert all([newname != f for f in folders])
	newname = os.path.abspath(newname)
	if os.path.exists(newname): 
			shutil.rmtree(newname)
	os.mkdir(newname)
	
	for fo in folders:
		flist = glob.glob(os.path.join(fo,"*.pdb"))
		if n!=None and len(flist)>n:
			flist = random.sample(flist, n)
		
		for fi in flist:
			
			flink = os.path.abspath(os.path.join(newname, os.path.basename(fi)))
			fsrc = os.path.abspath(fi)
			#print flink,fsrc
			if not os.path.exists(flink):
				if os.path.lexists(flink):
					os.remove(flink)
				os.symlink(fsrc, flink)
	return newname

def writeData2Pdb(pdb2data, labels, fname):
	df = open(fname,'w')
	df.write(",".join(labels)+"\n")
	for i in sorted(pdb2data):
		df.write("%s,%s\n"%(i,",".join([str(j) for j in pdb2data[i]])))
	df.close()

