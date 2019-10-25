import sys, math, os, argparse
import matplotlib.pyplot as p
from libprofsbm import get3Dhisto, findMinimumInHisto, wham, specificHeat, etotPerTemp, tempsFromDF, polyfitCv

import pandas as pd
import numpy as np

def run_anatraj(arguments):
		
	##################
	# HANDLE OPTIONS (-h for help)
	description = "%prog: apply WHAM to parallel tempering simulation data obtained from exTraj.py and plot 2D free energy surface."
	op = argparse.ArgumentParser(description=description)
	
	op.add_argument('csv_file', help="Input data (CSV) file, e.g. obtained by exTraj.py.")
	op.add_argument('indexA', help="Index/label of data column A")
	op.add_argument('indexB', help="Index/label of data column B")
	op.add_argument('-b','--begincycle', type=int, default=0, help="Ignore data before this MC cycle.")
	op.add_argument('-e','--endcycle', type=int, default=-1, help="Ignore data after this MC cycle.")
	op.add_argument('-s','--showPlot', action="store_true", help="Plot Cv curve and display to screen.")
	op.add_argument('-t','--targetTemp', type=float, help="Set target temperature other than Tm (default).")
	op.add_argument('-n','--normalizeBy', help="Specify normalization factors as N1,N2 so that the first (second) axis is divided by N1 (N2).")
	op.add_argument('-p','--polyfitDegree', type=int, default=0, help="Perform polynomial fit with given degree. ")
	op.add_argument('--min', help="Set axes minima as: minA,minB")
	op.add_argument('--max', help="Set axes maxima as: maxA,maxB")

	args = op.parse_args(arguments)
	
	df = pd.read_csv(args.csv_file)
	nrows = len(df.index)
	df = df.dropna()
	nnan = nrows - len(df.index)
	print "Dropped %i rows containing 'NaN'"%(nnan)
	
	assert all( [len(df.ix[:,0])==len(df.ix[:,i]) for i in xrange(len(df.columns))]  )
	
	A = df[args.indexA]
	B = df[args.indexB]
	E = df["Etot"]
	
	assert len(A)==len(B)==len(E)
	
	
	if args.normalizeBy:
		tmp = args.normalizeBy.split(",")
		assert len(tmp) == 2
		
		if tmp[0] == "": n1 = max(df[args.indexA])
		elif tmp[0] == "-": n1 = -max(df[args.indexA])
		else: n1 = float(tmp[0])
		
		if tmp[1] == "": n2 = max(df[args.indexB])
		elif tmp[1] == "-": n2 = -max(df[args.indexB])
		else: n2 = float(tmp[1])
		
		args.indexA = "%s/%f"%(str(args.indexA), n1)
		args.indexB = "%s/%f"%(str(args.indexB), n2)
		
		df[args.indexA] = A / n1
		df[args.indexB] = B / n2
		
		A = df[args.indexA]
		B = df[args.indexB]
	
	temps, tfreqs = tempsFromDF(df)
	print pd.DataFrame({"temps":temps, "freqs":tfreqs})
	
	if args.min: 
		QAmin, QBmin = map(float,args.min.split(","))
	else:
		QAmin = min(A)
		QBmin = min(B)
	
	if args.max: 
		QAmax, QBmax = map(float,args.max.split(","))
	else:
		QAmax = max(A)
		QBmax = max(B)
	
	Etotmin = min(E)
	Etotmax = max(E)
	
	QAlabel = args.indexA
	
	QBlabel = args.indexB
	
	print df.describe()
	#df = df.dropna(how='any')
	
	#print df.describe()
	
	print "NA count:", df["Etot"].hasnans()
	
	#SETTINGS
	nEtotbins = 40
	nQAbins = 101#51#101
	nQBbins = 102#52#102
	
	ncontours = nEtotbins#103#24
	
	qa_bisect = 0.5
	qb_bisect = 0.5
	
	maxFreeEnergy = 6.0##6.0#12.0#6.0#12.0#
	
	showPlots = False
	
	plotStyle = "default"; assert plotStyle in ["bare","default"]
	BGcolor = (1,1,1,0)
	if plotStyle=="bare": # just a square, no axes, no labels...
		showColorbar = False
		showCticks = True
		showCticklabels = True
		showClabel = False
		showXlabels = False
		showYlabels = False
		showXticks = True
		showYticks = True
		showXticklabels = False
		showYticklabels = False
		showAxes = False
	
		drawMainTitle = False 
		drawSubTitle = False
		drawFminX = False
		drawGuides = False
		drawFvalues = False
		drawDiagonal = False
		aspect = 'equal'
		
	elif plotStyle == "default":
		showColorbar = True
		showCticks = True
		showCticklabels = True
		showClabel = True
		showXlabels = True
		showYlabels = True
		showXticks = True
		showYticks = True
		showXticklabels = True
		showYticklabels = True
		showAxes = True
	
		drawMainTitle = True 
		drawSubTitle = True
		drawFminX = True
		drawGuides = False
		drawFvalues = True
		drawDiagonal = False
		aspect = 'equal'
	
	#Cv options
	
	minTmCut = 0.5
	maxTmCut = 1.0
	
	max_wham_iterations = 200
	useDecimalPackage = False # if getting OverflowError; is slow
	verbose = False
	figdpi = 400
	exportFmatrix = True
	
	
	#############################
	# INITIALIZATION
	
	qa_bisect = qa_bisect*QAmax
	qb_bisect = qb_bisect*QBmax

	QAbins = [QAmin + i*(QAmax-QAmin)/nQAbins for i in xrange(nQAbins ) ]
	QBbins = [QBmin + i*(QBmax-QBmin)/nQBbins for i in xrange(nQBbins ) ]
	
	print min(QAbins), max(QAbins)
	
	QAbinsize = QAbins[1]-QAbins[0]
	QBbinsize = QBbins[1]-QBbins[0]
	Ebinsize = (Etotmax-Etotmin)/nEtotbins
	Etotbins = [Etotmin + i*Ebinsize for i in xrange(nEtotbins ) ]
	
	QAmidpoints = [QAbins[i] + (QAbinsize/2.0) for i in xrange(nQAbins )  ]
	QBmidpoints = [QBbins[i] + (QBbinsize/2.0) for i in xrange(nQBbins )  ]
	Emidpoints = [Etotbins[i] + (Ebinsize/2.0) for i in xrange(nEtotbins )  ]
	
	assert QAmax-QAmin - len(QAbins)*QAbinsize < 10e-10, "%s"%(str(QAmax-QAmin-len(QAbins)*QAbinsize))
	assert QBmax-QBmin - len(QBbins)*QBbinsize < 10e-10
	assert Etotmax-Etotmin - len(Etotbins)*Ebinsize < 10e-10
	
	print "A:", QAmin,QAmax, len(QAbins), QAbinsize
	print "B:", QBmin,QBmax, len(QBbins), QBbinsize
	print "E:", Etotmin, Etotmax, len(Etotbins), Ebinsize
	
	if verbose:
		print "Etot min,max:",Etotmin,Etotmax
		print "Ebinsize:",Ebinsize
		print "Etotbins:",Etotbins
		print "Emidpoints:",Emidpoints

	print "Getting histogram..."
	# HISTOGRAMS
	H, histokeys = get3Dhisto(df, args.indexA, args.indexB, QAbins, QBbins, Etotbins, temps, mintime=args.begincycle, maxtime=args.endcycle)
	
	
	if args.targetTemp == None:
		print "Specific heat curve..."
		cvs, Tmelt, maxcv = specificHeat(temps, etotPerTemp(df, temps), maxTmCut, minTmCut)
		print cvs, Tmelt, maxcv
		if args.polyfitDegree >0:
			print "Specific heat (polynomial fit)..."
			fitT, fitcvs, Tmelt, fitMaxCv  = polyfitCv(temps, cvs, args.polyfitDegree, maxTmCut, minTmCut)
		
			print fitT, fitcvs, Tmelt, fitMaxCv
		args.targetTemp = Tmelt
	
	print "Applying WHAM..."
	P = wham(H, args.targetTemp, tfreqs, Emidpoints, max_wham_iterations, histokeys, temps, useDecimalPackage=useDecimalPackage, verbose=False, wham_cut=1e-3)
	
	# plot  2D surface wham
	
	fig, ax = p.Figure(figsize=(15, 15),dpi=800)
	bg = fig.patch
	bg.set_facecolor(BGcolor)
	title = os.path.basename(args.csv_file).replace("_rt.txt","")+";T%i-%i"%(args.begincycle,args.endcycle)
	
	if drawMainTitle:
		p.figtext(0.5, 0.965, title,ha='center', color='black', weight='normal', size='medium')

	
	ax.set_axis_bgcolor(BGcolor)
	ax.set_aspect(aspect)
	if not showAxes:
		ax.set_axis_off()
	
	if drawSubTitle:
		ax.set_title("Tmodel = %s "%(str(round(args.targetTemp,4) ) ) )
	p.xlabel(QBlabel)
	

	p.ylabel(QAlabel)
	
	
	#print P_unb_norm[t]
	summedEbins = [ [ sum([ P[(q1,q2,e)] if (q1,q2,e) in P else 0.0 for e in xrange(nEtotbins)] )  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
	#print summedEbins
	
	args.targetTemp
	
	minimum = min([ min([ -args.targetTemp * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else 1e900 for q2 in xrange(nQBbins) ]) for q1 in xrange(nQAbins)])
	maximum = max([ max([ -args.targetTemp * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else 1e-900 for q2 in xrange(nQBbins) ]) for q1 in xrange(nQAbins)])
	
	print minimum,maximum
	if maxFreeEnergy == None:
		dat = [ [ -args.targetTemp * math.log(summedEbins[q1][q2])-minimum if summedEbins[q1][q2] != 0.0 else maximum-minimum  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
	else:
		dat = [ [ -args.targetTemp * math.log(summedEbins[q1][q2])-minimum if summedEbins[q1][q2] != 0.0 else maxFreeEnergy  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
	#print dat
	
	#datB = [ [ -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else 0.0  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
	
	if exportFmatrix:
		Fname = os.path.basename(args.csv_file).replace("_rt.txt","_F_%s_%s.dat"%(QAlabel,QBlabel))
		print Fname
		outF = open(Fname,"w")
		
		nanDat = dat
		#nanDat = [ [ -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else None  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
		
		for i in nanDat:
			outF.write(",".join(map(str,i))+"\n")
		outF.close()
	
	minQA, minbinA_QA, minbinB_QA = findMinimumInHisto(dat, QAbins, QBbins, qa_bisect, max(1.0,QAmin), 0.0, qb_bisect)
	minQB, minbinA_QB, minbinB_QB = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, qb_bisect, max(1.0,QBmin))
	minUnf, minbinA_Unf, minbinB_Unf = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, 0.0, qb_bisect)
	
	print "Temperature",args.targetTemp
	print "Minimum at A:",minQA, minbinA_QA, minbinB_QA
	print "Minimum at B:",minQB, minbinA_QB, minbinB_QB
	print "Minimum at U:",minUnf, minbinA_Unf, minbinB_Unf
	print "FA-FB =",minQA-minQB
	print "FA-U =",minQA-minUnf
	print "FB-U =",minQB-minUnf
	print
	#print countPerQuadrantPerT(t,[datB], QAbins, QBbins, qa_bisect, qb_bisect,verbose=False)[0], "T%i ="%t,finalTemps[t]
	
	if drawDiagonal:
		p.plot([QBmin,QBmax],[QAmin,QAmax],marker="",color="0.25",linestyle=":", zorder=1)
	
	if maxFreeEnergy != None:
		levels = [0.0+i*(maxFreeEnergy/ncontours) for i in xrange(ncontours)]
		p.contourf(QBbins, QAbins, dat, ncontours,levels=levels, extend="max",figure=fig, zorder=2)
	else:
		p.contourf(QBbins, QAbins, dat , ncontours,figure=fig, zorder=2)
	
	if drawFvalues:
		if minbinB_QA!=None and minbinA_QA!=None: 
			p.text(QBbins[minbinB_QA]+0.5*QBbinsize, QAbins[minbinA_QA]+0.5*QAbinsize, str(round(minQA,3)), fontsize=12)
		if minbinB_QB!=None and minbinA_QB!=None:
			p.text(QBbins[minbinB_QB]+0.5*QBbinsize, QAbins[minbinA_QB]+0.5*QAbinsize, str(round(minQB,3)), fontsize=12)
		if minbinB_Unf!=None and minbinA_Unf!=None:
			p.text(QBbins[minbinB_Unf]+0.5*QBbinsize, QAbins[minbinA_Unf]+0.5*QAbinsize, str(round(minUnf,3)), fontsize=12)
	
	if drawFminX:
		if min(minQA,minQB,minUnf)==minQA and minbinB_QA!=None and minbinA_QA!=None:
			p.plot([QBbins[minbinB_QA]+0.5*QBbinsize],[QAbins[minbinA_QA]+0.5*QAbinsize],marker="x",color="0.5",markersize=10, markeredgewidth=3)
		elif min(minQA,minQB,minUnf)==minQB and minbinB_QB!=None and minbinA_QB!=None:
			p.plot([QBbins[minbinB_QB]+0.5*QBbinsize],[QAbins[minbinA_QB]+0.5*QAbinsize],marker="x",color="0.5",markersize=10, markeredgewidth=3)
		elif min(minQA,minQB,minUnf)==minUnf and minbinB_Unf!=None and minbinA_Unf!=None:
			p.plot([QBbins[minbinB_Unf]+0.5*QBbinsize],[QAbins[minbinA_Unf]+0.5*QAbinsize],marker="x",color="0.5",markersize=10, markeredgewidth=3)
	
	if drawGuides:
		p.plot(QBbins, [qa_bisect for b in QBbins],marker="",linewidth=0.25,color="k")
		p.plot([qb_bisect for b in QAbins],QAbins,marker="",linewidth=0.25,color="k")
	
	if showColorbar:
		if maxFreeEnergy != None:
			cbar = p.colorbar(ticks=[i for i in xrange(int(math.floor(maxFreeEnergy)))])
		else:
			cbar = p.colorbar()
	
		cbar.set_label("-T*ln(P)")
	
		if not showCticks:
			cbar.set_ticks([])
		if not showCticklabels:
			cbar.set_ticklabels([])
		if not showClabel:
			cbar.set_label("")
	
	if not showXlabels:
		ax.set_xlabel("")
	if not showYlabels:
		ax.set_ylabel("")
	if not showXticks:
		ax.set_xticks([])
	if not showYticks:
		ax.set_yticks([])
	if not showXticklabels:
		ax.set_xticklabels([])
	if not showYticklabels:
		ax.set_yticklabels([])
	
	figname = os.path.basename(args.csv_file).replace("_rt.txt","_F_%s_%s.png"%(QAlabel,QBlabel))
	
	print figname
	
	
	if not showAxes:
		p.axis('off')
	
	p.savefig(figname, dpi=figdpi,facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight', pad_inches = 0)	
	if showPlots: p.show()
	return 

if __name__=="__main__":
	print run_anatraj(sys.argv[1:])
