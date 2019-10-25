import sys, math, os, argparse
from scipy.optimize import curve_fit
import matplotlib.pyplot as p
from libprofsbm import geomTemp, optTbyOscillator, optTbyCurve, cooperativity, gaussfunc, gaussfuncScalar,tempsFromDF, specificHeat, etotPerTemp, polyfitCv

import pandas as pd
import numpy as np

def nReplicas(df):
	return df.replica.nunique()

def tempPerReplica(df):
	temps = sorted(df.temperature.unique())
	n = nReplicas(df)
	#print temps
	tcpr = np.array([[0 for t in xrange(df.temperature.nunique())] for r in xrange(n)])
	#print tcpr.shape
	for r in xrange(n):
		for t in xrange(len(temps)):
			c = len(df[df.replica==r][df.temperature==temps[t]])
			#print r,t, c
			tcpr[r][t] = c
	
	return tcpr
	




def swapProb(temps):
	p = []
	for t in xrange(len(temps)):
		for i in xrange(len(temps)):
			f = math.exp((temps[t] - temps[i]) / temps[t])
			if t < len(temps) - 1 and temps[i] == temps[t + 1]:
				p.append(f)
	return p

def run_heat(arguments):
	
	##################
	# HANDLE OPTIONS (-h for help)
	description = "%prog: plot specific heat (Cv) curve for a multi-temperature simulation."
	op = argparse.ArgumentParser(description=description)
	
	op.add_argument('csv_file', help="Input data (CSV) file, e.g. obtained by exTraj.py.")
	op.add_argument('-b','--begincycle', type=int, default=0, help="Ignore data before this MC cycle.")
	op.add_argument('-e','--endcycle', type=int, default=-1, help="Ignore data after this MC cycle.")
	op.add_argument('-s','--showPlot', action="store_true", help="Plot Cv curve and display to screen.")
	op.add_argument('-p','--polyfitDegree', type=int, default=0, help="Perform polynomial fit with given degree. Default=%default")
	args = op.parse_args(arguments)
	
	df = pd.read_csv(args.csv_file)
	
	print df.describe()
	
	#Cv options
	minTmCut = 0.1
	maxTmCut = 1.0
	
	tgrid_adjustRange = False
	autoSuggestTrange = False
	autoTdelta = 20
	newTmin = 0.05
	newTmax = 0.05
	relative_to_Tm = True
	set_ntemp = 32
	baseGridGeometric = True
	optTmethod = "none" ; assert optTmethod in ["slope","curve","none"]
	
	vantHoff = False
	
	
	figlowdpi = 100
	
	#############################
	
	# INITIALIZATION
	
	temps, freqs = tempsFromDF(df)
	print pd.DataFrame({"temps":temps, "freqs":freqs})
	
	
	
	cvs, Tmelt, maxcv = specificHeat(temps, etotPerTemp(df, temps), maxTmCut=maxTmCut, minTmCut=minTmCut)
	print "Tm=%f"%Tmelt
	print "max(Cv)=%f"%maxcv
	
	if relative_to_Tm:
		newTmin = Tmelt-newTmin
		newTmax = Tmelt+newTmax
	
	sp = swapProb(temps)
	tempCountPerReplica = tempPerReplica(df)
	print tempCountPerReplica
# 	fig = p.figure(figsize=(7, 3))
# 	ax = fig.add_subplot(111)
# 	H, xedges, yedges = np.histogram2d(df["replica"],df["temperature"], bins=len(temps))
# 	ax.set_title('NonUniformImage: interpolated')
# 	import matplotlib as mpl
# 	im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
# 	xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
# 	ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
# 	im.set_data(xcenters, ycenters, H)
# 	ax.images.append(im)
# 	ax.set_xlim(xedges[0], xedges[-1])
# 	ax.set_ylim(yedges[0], yedges[-1])
# 	ax.set_aspect('auto')
# 	
# 	p.show()
# 	
# 	sys.exit()
	print "Swap prob:",sp
 	print "\nTemperatures visited per replica:"
 	print  "T:\t" + "\t".join([str(round(i,3)) for i in temps])
 	for r in xrange(len(tempCountPerReplica)):
 		print "R%i:\t"%r + "\t".join([str(tempCountPerReplica[r][t]) if t in tempCountPerReplica[r] else "0"  for t in reversed(xrange(len(temps)))])
 	print "\nCv(T):\t"  + "\t".join([str(round(i,3)) for i in cvs])
	
	
	
	#print "\nTm=%f,maxcv=%f\n"%(Tmelt,maxcv)
	
	if autoSuggestTrange:
		newTmin = Tmelt-autoTdelta
		newTmax = Tmelt+autoTdelta
	
	if tgrid_adjustRange or autoSuggestTrange:
		
		if newTmax==None:
			newTmax = Tmelt+(Tmelt-newTmin)
		
		if baseGridGeometric:
			Ttmp = geomTemp(newTmin,newTmax,set_ntemp)
		else:
			Ttmp = [newTmin + (i*(newTmax-newTmin)/(set_ntemp-1)) for i in xrange(set_ntemp)]
	else:
		if baseGridGeometric:
			Tnew = geomTemp(min(temps),max(temps),len(temps))
		else:
			Tnew = [min(temps) + (i*(max(temps)-min(temps))/len(temps)) for i in xrange(len(temps))]
		Ttmp = None
	
	if optTmethod == "slope":
		Tnew, cvpeaks = optTbyOscillator(temps,cvs, rt=5.0, adjustRange=Ttmp)
	elif optTmethod == "curve":
		if args.polyfitDegree>0:
			Tnew, cvpeaks, newC = optTbyCurve(temps,cvs, rt=5.0, adjustRange=Ttmp,polyfit=args.polyfitDegree)
		else:
			Tnew, cvpeaks, newC = optTbyCurve(temps,cvs, rt=5.0, adjustRange=Ttmp)
	elif optTmethod == "none":
		if Ttmp != None: Tnew = Ttmp
	
	#sys.exit()
	if vantHoff:
		#TODO: get begining and end of peak from 2nd derivative of polyfit curve
		left_end_index = 32
		right_begin_index = 120
		coop, baselineCorrectedCv, left_linCv, right_linCv = cooperativity(temps,cvs,Tmelt,maxcv, left_end_index=left_end_index, right_begin_index=right_begin_index)
	
		print np.array(temps)
		print np.array(baselineCorrectedCv)
		gauss_popt, gauss_pcov = curve_fit( gaussfunc, np.array(temps), np.array(baselineCorrectedCv),p0=np.array([maxcv,1]))
	
		print "gauss_popt",gauss_popt
		print "gauss_pcov",gauss_pcov
	
		print left_end_index,left_linCv
		print right_begin_index, right_linCv
		print "van't Hoff calorimetric cooperativity: ", coop
	
	
	
 	print  "Told:\t" + "\t".join([str(round(i,3)) for i in temps])
	print  "Tnew:\t" + "\t".join([str(round(i,3)) for i in Tnew])

	print"[begin COPY&PASTE for T.dat]"
	print "#temperature"
	for i in reversed(Tnew):
		print i
	print"[end COPY&PASTE for T.dat]"
	
	p.clf()
	fig, ax1 = p.subplots()
	fig.set_size_inches(18.5, 10.5)		
	
	ax1.set_title("%s; MCcyc:%sM to %sM\nTm=%s"%(os.path.basename(args.csv_file), str(round(args.begincycle/1000000.0,3)),str(round(args.endcycle/1000000.0,3)), str(round(Tmelt,4))))
	ax1.set_ylabel("Cv(T)")
	ax1.set_xlabel("T(model)")
	
	ax1.plot(temps,cvs,marker="x",label="Cv",linewidth=2.0)
	
	print "Tm=%f; Cv(Tm)=%f"%(Tmelt, maxcv)
	
	if args.polyfitDegree >0:
		fitT, fitcvs, fitTmelt, fitMaxCv  = polyfitCv(temps, cvs, args.polyfitDegree, maxTmCut, minTmCut)
		ax1.plot(fitT,fitcvs,marker="",label="fitted Cv (deg %i)"%args.polyfitDegree,linewidth=2.0)
		#ax2 = ax1.twinx()
		#ax2.plot(finegrid,fineDCv,marker="",label="fitted Cv (deg %i)"%polyfitDegree,linewidth=1.0)
		#ax2.set_label("2nd derivative")
		print "Polynomial fit to Cv curve, degree=",args.polyfitDegree
		print "fitted Tm, max_Cv:",fitTmelt,fitMaxCv
		ax1.set_title(ax1.get_title()+" (polyfit(%i) Tm=%s)"%(args.polyfitDegree, str(round(fitTmelt,4))))
	print "Tm forced to be between %f and %f"%(minTmCut, maxTmCut)
	
	if vantHoff:
		ax1.plot(temps,baselineCorrectedCv,marker="+",label="Cv (for coop.)",linewidth=2.0)
		ax1.plot(temps,left_linCv,marker=".",label="left baseline",linewidth=2.0)
		ax1.plot(temps,right_linCv,marker=".",label="right baseline",linewidth=2.0)
		ax1.plot([temps[left_end_index],temps[left_end_index]],[0,maxcv],label="T0",marker="", color="r",linewidth=1,linestyle=":")
		ax1.plot([temps[right_begin_index],temps[right_begin_index]],[0,maxcv],label="T1",marker="", color="r",linewidth=1,linestyle=":")
	
		g = [gaussfuncScalar(c, gauss_popt[0], gauss_popt[1],m=Tmelt) for c in baselineCorrectedCv]
		print g, "g"
		ax1.plot(temps,g,marker="",label=" Gauss of Cv (for coop.)",linewidth=2.0)
	
	if optTmethod == "curve":
		ax1.plot(Tnew,newC,marker="o",label="Cv",linewidth=2.0)
	
	if tgrid_adjustRange:
		ax1.plot(temps,[-200 for i in xrange(len(temps))],marker="o",label="Told")
		ax1.plot(Tnew,[-100 for i in xrange(len(Tnew))],marker="v",label="Tnew")
		ax1.set_ylim(-300)
	
	
	ax1.legend(loc=0)
	p.savefig(os.path.basename(args.csv_file).replace(".csv","_Cv.png"), dpi=figlowdpi)
	
	#write Cv data
	cvout = open(os.path.basename(args.csv_file).replace(".csv","_Cv.csv"),"w")
	cvout.write(",".join(["T","Cv"])+"\n")
	for i in xrange(len(temps)):
		cvout.write(",".join(map(str,[temps[i],cvs[i]]))+"\n")
	cvout.close()
	
	if args.polyfitDegree >0:
		#write fitted Cv data
		cvout = open(os.path.basename(args.csv_file).replace(".csv","_CvFit.csv"),"w")
		cvout.write(",".join(["T","fitted_Cv"])+"\n")
		for i in xrange(len(fitT)):
			cvout.write(",".join(map(str,[fitT[i],fitcvs[i]]))+"\n")
		cvout.close()
	#sys.exit()
	if args.showPlot: p.show()
	return cvout.name


if __name__=="__main__":
	print run_heat(sys.argv[1:])
