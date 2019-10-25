import sys, random, math, numpy, os, operator, shutil, optparse
from scipy.optimize import curve_fit
import matplotlib.pyplot as p
import matplotlib.patches as mpatches
from libprofsbm import checkHeader, getTemps, autoCorrelation, getCv, doConstTempStats, parseDataFromSingleFile, getHistograms, KtoPRFbeta, findMinimumInHisto, wham, swapProb, getDiffusivity, geomTemp, optTbyOscillator, optTbyCurve, cooperativity, gaussfunc, gaussfuncScalar, isInState, findMinimumInHisto1D, getBinIndex, getCoordPDB, isPointInRectangles

#import pandas as pd

def run_anatraj(opts, arguments):
	assert len(arguments) >= 1
	datafile = arguments[0]
	#df = pd.read_csv(datafile)
	#print df.keys()
	if len(arguments) >= 6:
		
		nativenessQ = arguments[3]
	
		indexA = int(arguments[4])
		indexB = int(arguments[5])
	else:
		nativenessQ = None
		indexA = None
		indexB = None
	
	if len(arguments) >= 7:
		EAmin = float(arguments[6])
		EBmin = float(arguments[7])
	else:
		EAmin = None
		EBmin = None
	
	headers, datmins, datmaxs, datavg, datstd = checkHeader(datafile)
	if headers != []:
		labeldict = {}
		print "i","term", "min","max","avg","std"
		for i in xrange(len(headers)):
			labeldict[i] = headers[i]
			print i, labeldict[i], datmins[i], datmaxs[i], datavg[i], datstd[i]
		
	else:
		labeldict={
		0	:"Time in MC cycles",
		1	:"temperature",
		2	:"Etot",
		3	:"ExVol",
		4	:"LocExVol",
		5	:"Bias",
		6	:"TorsionTerm",
		7	:"HBMM",
		8	:"HBMS",
		9	:"Hydrophobicity",
		10:"ChargedSCInteraction",
		11:"DistanceRestraints_superGA",
		12:"DistanceRestraints_superGB",
		13:"rg",
		14:"rmsd_2LHC",
		15:"rmsd_2LHD",
		16:"Q_2LHC",
		17:"Q_2LHD",
		}
		for i in xrange(len(headers)):
			print i, labeldict[i]
	
	if nativenessQ == None or indexA==None or indexB==None: sys.exit()
	
	#SETTINGS
	
	minsteps = 3000000 #500000#0#3000000
	maxsteps = 10000000 #6500000#10000000
	stepsize = 1
	
	ncontours = 103#24
	
	nEtotbins = 40
	nQAbins = 101#51#101
	nQBbins = 102#52#102
	
	Tmix = None
	weightA = 1.0
	weightGO = 1.0
		
	if Tmix != None: Tmix = Tmix * weightGO
	
	print_free_energy = False
	print_histogram = False
	
	qa_bisect = 0.5
	qb_bisect = 0.5
	
	maxFreeEnergy = 6.0##6.0#12.0#6.0#12.0#
	
	showPlots = False
	
	plotStyle = "default"; assert plotStyle in ["bare","default"]
	
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
	
	Tunit = "Profasi"; assert Tunit in ["Kelvin","Profasi"]
	#Cv options
	
	minTmCut = 0.5
	maxTmCut = 1.0
	
	max_wham_iterations = 200
	whamMinP = 1e-99
	
	conCheckTrajChunks = 12
	checkConvCumulative = True
	
	tgrid_adjustRange = False
	autoSuggestTrange = False
	autoTdelta = 20
	newTmin = 0.05
	newTmax = 0.05
	relative_to_Tm = True
	set_ntemp = 32
	baseGridGeometric = True
	optTmethod = "none" ; assert optTmethod in ["slope","curve","none"]
	
	useDecimalPackage = False # if getting OverflowError; is slow
	
	verbose = False
	
	vantHoff = False
	
	polyfitDegree = 3   # None to turn off Cv curve polynomial fitting, otherwise polynomial degree
	
	figdpi = 400
	figlowdpi = 100
	
	exportFmatrix = True
	
	checkAutoCorrelation = False
	checkConvergence = False
	trackStateTransitions = False
	analyseTempDiffusion = True
	constTempStats = False
	freeEnergy1D = True
	freeEnergy2D = True  
	lundVsGo = True
	
	############################# CLUSTERS:
	markStructuresFromFile = False
	
	if markStructuresFromFile: 
		ncontours = 2
		
	feature = "superposition"; assert feature in ["superdihedralcontacts","contacts","dihedrals","superposition","rmsd"]
	drawClusterMembers = True
	superpos_refpdbs = []
	superpos_atypes = ["CA", "CB"]
	if feature=="contacts":
		centroidDistCutoff = 5.75
		#minClusterSizeFraction = 0.01
		minClusterSizeFraction = 0.0
		sizeboost=1.0
		connectStrength = 1.0
	elif feature=="dihedrals":
		centroidDistCutoff = 0.575
		#minClusterSizeFraction = 0.005
		minClusterSizeFraction = 0.0
		sizeboost = 1.0
		connectStrength = 3.0
	elif feature=="superposition":
		superpos_refpdbs = ["2LHC_m1.pdb", "2LHD_m1.pdb"]
		centroidDistCutoff = 0.575
		#minClusterSizeFraction = 0.005
		minClusterSizeFraction = 0.0
		sizeboost = 1.0
		connectStrength = 1.5
	elif feature == "superdihedralcontacts":
		superpos_refpdbs = ["2LHC_m1.pdb", "2LHD_m1.pdb"]
		centroidDistCutoff = 0.575
		#minClusterSizeFraction = 0.005
		minClusterSizeFraction = 0.0
		sizeboost = 1.0
		connectStrength = 3.0
		
	drawCentroidCircles = True
	drawCentroidLabels = False
	drawCentroidConnections = True
	drawCentroidLegend = False
	drawCentroidRandomColors = False
	centroidRandomColorsSeed = 123 # any hashable or None
	random.seed(centroidRandomColorsSeed)
	uniformCentroidColor = (0.5,0.5,0.5,1.0)
	BGcolor = (1,1,1,0)
	highlightclusters = []#[43,46,23,0,25,12,22,29,13]
	cluster_highlight_color = (1,1,0,0.5)
	highlightNtopClusters = 3
	minCentroidConn = 2
	
	minClusterIndex = 0 # start cluster index with ...
	
	onlyDrawMembersInRectangles = True
	analyseRectangles = True
	#TS GA
	#qa_cluster_range = (0.66, 0.74)  # only draw cluster members in this range of coordinates
	#qb_cluster_range = (0.12, 0.22)

	#TS GB1
	#qa_cluster_range = (0.30, 0.55)  # only draw cluster members in this range of coordinates
	#qb_cluster_range = (0.35, 0.43)

	#TS GB2
	#qa_cluster_range = (0.28, 0.4)  # only draw cluster members in this range of coordinates
	#qb_cluster_range = (0.58, 0.66)
	
	drawRectangles = [
		(0.12, 0.22, 0.66, 0.74, cluster_highlight_color, "TS-GA"),
		(0.35, 0.43, 0.30, 0.55, cluster_highlight_color, "TS1-GB"),
		(0.58, 0.66, 0.28, 0.4, cluster_highlight_color, "TS2-GB")
	]
	
	exportPDBsInRange = False
		
	#centroids_fname = "/Users/sikosek/protein/MC_results/GB30_sAB3.75_A48.8_B51_fineTm_ALL_all_N20000_kmeans50.txt"
	#centroids_fname = "/Users/sikosek/protein/MC_results/ProtL_sAB3.75_A48.8_B51_fineTm_ALL_all_N20000_kmeans100_contacts.txt"
	#GBwt_sAB3.75_A48.8_B51_fineTm_ALL_all_N20000_kmeans50_dihedrals.txt"
	#"Gx98_N2500m_kmeans20_dihedrals.txt"
	#centroids_fname = "/Users/sikosek/protein/MC_results/GB98_A48.8_opt_ALL_all_N20000_kmeans100_contacts.txt"
	centroids_fname = "/Users/sikosek/protein/MC_results/Gx98_sAB0_A48.8_B51_fineTm_ALL_all_N20000_merged_kmeans50_superposition.txt"
	#"/Users/sikosek/protein/MC_results/Gx98_k100means_centroids_maxcluster_C3_is_15.txt"
	#
	#Gx98_sAB3.75_A48.8_B51_fineTm_ALL_all_N2500_merged_is283_centroids.txt" 	
	#"/Users/sikosek/protein/MC_results/GB98_sAB3.75_A48.8_B51_fineTm_ALL_unfolded_N1000/centroids.txt"
	#"/Users/sikosek/protein/MC_results/GB98_sAB3.75_A48.8_B51_fineTm_ALL_all_N2500/centroids.txt"
	#"/Users/sikosek/protein/MC_results/GA98_sAB3.75_A48.8_B51_fineTm_ALL_all_N5000_maxcluster_centroids.txt"
	# "/Users/sikosek/protein/MC_results/GA98_sAB3.75_A48.8_B51_fineTm_ALL_all_N5000_maxcluster_3djury_top20.txt" #
	#pdb2data_fname =  "/Users/sikosek/protein/MC_results/GB30_sAB3.75_A48.8_B51_fineTm_ALL_all_N20000_pdb2data.dat"
	#pdb2data_fname =  "/Users/sikosek/protein/MC_results/GB98_A48.8_opt_ALL_all_N20000_pdb2data.dat"
	pdb2data_fname =  "/Users/sikosek/protein/MC_results/Gx98_sAB0_A48.8_B51_fineTm_ALL_all_N20000_merged_pdb2data.dat"
	#"/Users/sikosek/protein/MC_results/GB98_sAB3.75_A48.8_B51_fineTm_ALL_unfolded_N1000_pdb2data.dat"
	#"/Users/sikosek/protein/MC_results/GB98_sAB3.75_A48.8_B51_fineTm_ALL_all_N2500_pdb2data.dat"
	#pdb2data_fname =  "/Users/sikosek/protein/MC_results/ProtL_sAB3.75_A48.8_B51_fineTm_ALL_all_N20000_pdb2data.dat"
	#GBwt_sAB3.75_A48.8_B51_fineTm_ALL_all_N20000_pdb2data.dat"
	#Gx98_sAB3.75_A48.8_B51_fineTm_ALL_all_N2500_merged_pdb2data.dat"
	
	
	#############################
	
	# INITIALIZATION
	
	temps = getTemps(tempfile)
	if Tunit=="Kelvin":
		tkel = [temps[t][1] for t in reversed(xrange(len(temps)))]
	else:
		tkel = [temps[t][0] for t in reversed(xrange(len(temps)))]
	
	finalTemps = []#tkel#[]#[250]#[788,716,537,250]#range(333)# if empty: auto-generate
	
	if nativenessQ in ["True","true","1"]:
		if EAmin == None or EBmin == None:
			QAmin = datmins[indexA]
			QAmax = datmaxs[indexA]
			
			QBmin = datmins[indexB]
			QBmax = datmaxs[indexB]
		else:
			if EAmin > 0:
				QAmin = datmins[indexA]#QAmin = 0.0
				QAmax = EAmin +0.1#1.1
			else:
				QAmin = EAmin
				QAmax = datmaxs[indexA]#QAmax = 2.0
		
			if EBmin > 0:
				QBmin = datmins[indexB]#QBmin = 0.0
				QBmax = EBmin +0.1#1.1
			else:
				QBmin = EBmin
				QBmax = datmaxs[indexB]#QBmax = 2.0
			
		QAlabel = labeldict[indexA]# "RMSD_A"#"E_Hydrophobicity"#"RMSD A"#"QA"
		
		QBlabel = labeldict[indexB]#"RMSD_B"#"RMSD B"#"RMSD B"#
	else:
		if EAmin == None or EBmin == None:
			EAmin = datmins[indexA]
			EBmin = datmins[indexB]
		QAmin = 0.0#-0.1
		QAmax = 1.0#1.1
		QAlabel = "-EA/%s"%str(EAmin)
		QBmin = 0.0#-0.1
		QBmax = 1.0#1.1
		QBlabel = "-EB/%s"%str(EBmin)
	
	qa_bisect = qa_bisect*QAmax
	qb_bisect = qb_bisect*QBmax
	
	print QAmin,QAmax, QBmin,QBmax
	
	QAbins = [QAmin + i*(QAmax-QAmin)/nQAbins for i in xrange(nQAbins ) ]
	QBbins = [QBmin + i*(QBmax-QBmin)/nQBbins for i in xrange(nQBbins ) ]
	
	#print QAbins
	#print QBbins
	
	QAbinsize = QAbins[1]-QAbins[0]
	QBbinsize = QBbins[1]-QBbins[0]
	
	QAmidpoints = [QAbins[i] + (QAbinsize/2.0) for i in xrange(nQAbins )  ]
	QBmidpoints = [QBbins[i] + (QBbinsize/2.0) for i in xrange(nQBbins )  ]
	
	#print QAmidpoints
	#print QBmidpoints
	
	raw_data, etotPerTemp, egoPerTemp, nT0, nT1, tempCountPerReplica, Etotmin, Etotmax, maxtime, raw_data_per_replica, timeStep = parseDataFromSingleFile(datafile, temps, indexA, indexB, minsteps, maxsteps,stepsize, nativenessQ, EAmin, EBmin, constTempStats, Tmix=Tmix,weightA=weightA,weightGO=weightGO,lundVsGo=lundVsGo)
	
	Ebinsize = (Etotmax-Etotmin)/nEtotbins
	Etotbins = [Etotmin + i*Ebinsize for i in xrange(nEtotbins ) ]
	Emidpoints = [Etotbins[i] + (Ebinsize/2.0) for i in xrange(nEtotbins )  ]
	
	if verbose:
		print "Etot min,max:",Etotmin,Etotmax
		print "Ebinsize:",Ebinsize
		print "Etotbins:",Etotbins
		print "Emidpoints:",Emidpoints
	
	print [len(raw_data_per_replica[i]) for i in xrange(len(raw_data_per_replica))]
	
	if any([len(etotPerTemp[t])==0.0 for t in xrange(len(etotPerTemp))]):
		# data from constant temperature run. assuming replicas follow same order and index of temperatures
		
		for r in xrange(len(raw_data_per_replica)):
			etotPerTemp[r] = []
			
			for i in xrange(len(raw_data_per_replica[r])):
				E = float(raw_data_per_replica[r][i][2])
				etotPerTemp[r].append(E)
	
	#print [len(etotPerTemp[i]) for i in xrange(len(etotPerTemp))]
	
	
		
	if checkAutoCorrelation:
		if constTempStats:
			const_temps = {}
			for i in xrange(len(raw_data_per_replica)):
				const_temps[i] = (temps[0][0],temps[0][1],temps[0][2])
			const_etotPerTemp = [[raw_data_per_replica[i][j][2] for j in xrange(len(raw_data_per_replica[i]))]  for i in xrange(len(raw_data_per_replica))]
			print "MC sweeps per data step:",timeStep
			autoCorrelation(const_temps,const_etotPerTemp,datafile, maxtime, nSweepsPerData=timeStep, showPlots=showPlots, filename=os.path.basename(datafile).replace("_rt.txt","_ac.png"))
		else:
			autoCorrelation(temps,etotPerTemp,datafile, maxtime, nSweepsPerData=timeStep,showPlots=showPlots,filename=os.path.basename(datafile).replace("_rt.txt","_ac.png"))
		#sys.exit()
	
	if checkConvergence:
		cvs, Tmelt, Tmelti, maxcv, dCv  = getCv(temps,etotPerTemp,maxTmCut=maxTmCut,minTmCut=minTmCut, polyfitDegree=polyfitDegree, Tunit=Tunit)
		fixedminsteps = minsteps
		chunks = conCheckTrajChunks
		raw_data, etotPerTemp, egoPerTemp, nT0, nT1, tempCountPerReplica, Etotmin, Etotmax, maxtime, raw_data_per_replica, timeStep = parseDataFromSingleFile(datafile, temps, indexA, indexB, minsteps, maxsteps, stepsize, nativenessQ, EAmin, EBmin, constTempStats, Tmix=Tmix, weightA=weightA, weightGO=weightGO, lundVsGo=lundVsGo)
		timestep = (maxtime-minsteps)/float(chunks)
		fafb=[]
		fau=[]
		fbu=[]
		times=[]
		for i in xrange(chunks):
			if checkConvCumulative:
				minsteps = fixedminsteps
				maxsteps = minsteps+ ((i+1)*timestep)
				
			else:
				minsteps = i*timestep
				maxsteps = minsteps+timestep
			print minsteps,maxsteps,chunks,timestep
			raw_data, etotPerTemp, egoPerTemp, nT0, nT1, tempCountPerReplica, Etotmin, Etotmax, maxtime, raw_data_per_replica, timeStep = parseDataFromSingleFile(datafile, temps, indexA, indexB, minsteps, maxsteps, stepsize, nativenessQ, EAmin, EBmin, constTempStats, Tmix=Tmix, weightA=weightA, weightGO=weightGO, lundVsGo=lundVsGo)
			doConstTempStats(raw_data_per_replica)
			Ebinsize = (Etotmax-Etotmin)/nEtotbins
			Etotbins = [Etotmin + i*Ebinsize for i in xrange(nEtotbins ) ]
			Emidpoints = [Etotbins[i] + (Ebinsize/2.0) for i in xrange(nEtotbins )  ]
			histokeys, histokeysA, histokeysB, countPerTemp, histoPerTemp3Deff, totalCount, histoPerTempQAE, histoPerTempQBE, histoPerTempEtot,histoPerTempQ1,histoPerTempQ2 = getHistograms(raw_data, QAbins, QBbins, Etotbins, temps)
			finalTemps = [Tmelt]
			
			P_unb_norm = wham(finalTemps, histoPerTemp3Deff, countPerTemp, Emidpoints, max_wham_iterations, histokeys, temps, useDecimalPackage=useDecimalPackage, verbose=verbose, wham_cut=1e-3, Tunit=Tunit)
			for t in xrange(len(finalTemps)):
				summedEbins = [ [ sum([ P_unb_norm[t][(q1,q2,e)] if (q1,q2,e) in P_unb_norm[t] else 0.0 for e in xrange(nEtotbins)] )  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]

				maximum = max([ max([ -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else None  for q2 in xrange(nQBbins)]) for q1 in xrange(nQAbins)])

				dat = [ [ -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] > whamMinP else -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(whamMinP)  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
				#minQA, minbinA_QA, minbinB_QA = findMinimumInHisto(dat, QAbins, QBbins, qa_bisect, 1.0, 0.0, qb_bisect)
				#minQB, minbinA_QB, minbinB_QB = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, qb_bisect, 1.0)

				#print "Temperature",finalTemps[t]
				#print "Minimum at A:",minQA, minbinA_QA, minbinB_QA
				#print "Minimum at B:",minQB, minbinA_QB, minbinB_QB
				
				
				
				minQA, minbinA_QA, minbinB_QA = findMinimumInHisto(dat, QAbins, QBbins, qa_bisect, max(1.0,EAmin), 0.0, qb_bisect)
				minQB, minbinA_QB, minbinB_QB = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, qb_bisect, max(1.0,EBmin))
				minUnf, minbinA_Unf, minbinB_Unf = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, 0.0, qb_bisect)
				print "time(%i,%i)"%(minsteps,maxtime)
				print "Temperature",finalTemps[t]
				print "Minimum at A:",minQA, minbinA_QA, minbinB_QA
				print "Minimum at B:",minQB, minbinA_QB, minbinB_QB
				print "Minimum at U:",minUnf, minbinA_Unf, minbinB_Unf
				print "FA-FB =",minQA-minQB
				print "FA-U =",minQA-minUnf
				print "FB-U =",minQB-minUnf
				print
				times.append(maxtime)
				fafb.append(minQA-minQB)
				fau.append(minQA-minUnf)
				fbu.append(minQB-minUnf)
		convFile = open(datafile[:-4]+"_Fconvergence.csv",'w')
		print "last_time(MCcyc>=%i;cum?%s), FA(qa>%s)-FB(qb>%s), FA-U, FB-U"%(minsteps, str(checkConvCumulative), str(qa_bisect), str(qb_bisect))
		convFile.write("last_time(MCcyc>=%i;cum?%s), FA(qa>%s)-FB(qb>%s), FA-U, FB-U\n"%(minsteps, str(checkConvCumulative), str(qa_bisect), str(qb_bisect)))
		for f in xrange(len(fafb)):
			print times[f],fafb[f],fau[f],fbu[f]
			convFile.write(",".join(map(str,[times[f],fafb[f],fau[f],fbu[f]]))+"\n")
		convFile.close()
		#sys.exit()
	
	if analyseTempDiffusion:
		
		cvs, Tmelt, Tmelti, maxcv, dCv  = getCv(temps, etotPerTemp, maxTmCut=maxTmCut, minTmCut=minTmCut, Tunit=Tunit)
		
		if relative_to_Tm:
			newTmin = Tmelt-newTmin
			newTmax = Tmelt+newTmax
		
		sp = swapProb(tkel)
		
		print "Swap prob:",sp
		print "\nTemperatures visited per replica:"
		print  "T(K):\t" + "\t".join([str(round(i,3)) for i in tkel])
		for r in xrange(len(tempCountPerReplica)):
			print "R%i:\t"%r + "\t".join([str(tempCountPerReplica[r][t]) if t in tempCountPerReplica[r] else "0"  for t in reversed(xrange(len(temps)))])
		print "\nCv(T):\t"  + "\t".join([str(round(i,3)) for i in cvs])
		
		if False: # experimental feature: Duffusivity
			fT, prob, Diff = getDiffusivity(nT0,nT1,temps)
			psum = sum([ prob[i]*(tkel[i+1]-tkel[i]) for i in range(0,len(tkel)-1)])
			print "f(T):\t"   + "\t".join(str(round(i,3)) for i in fT)
			print "Ti+1/Ti:\t" + "\t".join([str(tkel[i+1] / tkel[i]) for i in range(0,len(tkel)-1)])
			print "P:",prob,psum
			print "Diff:",Diff
		
		print "\nTm=%f,maxcv=%f\n"%(Tmelt,maxcv)
		
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
				Tnew = geomTemp(min(tkel),max(tkel),len(tkel))
			else:
				Tnew = [min(tkel) + (i*(max(tkel)-min(tkel))/len(tkel)) for i in xrange(len(tkel))]
			Ttmp = None
		
		if optTmethod == "slope":
			Tnew, cvpeaks = optTbyOscillator(tkel,cvs, rt=5.0, adjustRange=Ttmp)
		elif optTmethod == "curve":
			if polyfitDegree:
				Tnew, cvpeaks, newC = optTbyCurve(tkel,cvs, rt=5.0, adjustRange=Ttmp,polyfit=polyfitDegree)
			else:
				Tnew, cvpeaks, newC = optTbyCurve(tkel,cvs, rt=5.0, adjustRange=Ttmp)
		elif optTmethod == "none":
			if Ttmp != None: Tnew = Ttmp
		
		#sys.exit()
		if vantHoff:
			#TODO: get begining and end of peak from 2nd derivative of polyfit curve
			left_end_index = 32
			right_begin_index = 120
			coop, baselineCorrectedCv, left_linCv, right_linCv = cooperativity(tkel,cvs,Tmelt,maxcv, left_end_index=left_end_index, right_begin_index=right_begin_index)
		
			print numpy.array(tkel)
			print numpy.array(baselineCorrectedCv)
			gauss_popt, gauss_pcov = curve_fit( gaussfunc, numpy.array(tkel), numpy.array(baselineCorrectedCv),p0=numpy.array([maxcv,1]))
		
			print "gauss_popt",gauss_popt
		
			print left_end_index,left_linCv
			print right_begin_index, right_linCv
			print "van't Hoff calorimetric cooperativity: ", coop
		
		if False:
			a = maxcv/1.0
			b = 0.8
			c = 10.0
			print "a=%d,b=%d"%(a,b)
	
			Tnew = [tkel[0], Tmelt, tkel[-1]]
	
			print "new temperatures below Tm:"
			print "abs(dC),  k =  b * math.sqrt((abs(dC)) / (a + abs(dC)))"
			Ttmp = Tmelt
			for i in range(Tmelti,1,-1):
				dC = (cvs[i]-cvs[i-1])
				#k = (abs(dC)) / (a + abs(dC))
				#k =  b * math.sqrt((abs(dC)) / (a + abs(dC)))
				k =  b * math.pow((abs(dC)) / (a + abs(dC)), 1.0/c)
				print abs(dC), k
				Ttmp = tkel[i-1] + (Ttmp-tkel[i-1])* (k)
				Tnew.append(Ttmp)
		
			print
			print "new temperatures above Tm:"
			print "abs(dC),  k =  b * math.sqrt((abs(dC)) / (a + abs(dC)))"
			Ttmp = Tmelt
			for i in range(Tmelti,len(tkel)-2,1):
				dC = (cvs[i+1]-cvs[i])
				#k = (abs(dC)) / (a + abs(dC))
				#k =  b * math.sqrt((abs(dC)) / (a + abs(dC)))
				k =  b * math.pow((abs(dC)) / (a + abs(dC)), 1.0/c)
				#k=0.75
				print abs(dC), k
				Ttmp = Ttmp + (tkel[i+1]-Ttmp)*(1.0 - k)
				Tnew.append(Ttmp)
	
			Tnew = sorted(Tnew)
			print 
			print "Swap prob Tnew:"
			sp = swapProb(Tnew)
	
			assert len(Tnew) == len(tkel)
		
	 	print  "Told:\t" + "\t".join([str(round(i,3)) for i in tkel])
		print  "Tnew:\t" + "\t".join([str(round(i,3)) for i in Tnew])
	
		print"[begin COPY&PASTE for T.dat]"
		if Tunit=="Kelvin":
			print "#temperature Kelvin"
		else:
			print "#temperature"
		for i in reversed(Tnew):
			print i
		print"[end COPY&PASTE for T.dat]"
	
		fig, ax1 = p.subplots()
		fig.set_size_inches(18.5, 10.5)		
		
		ax1.set_title("%s; MCcyc:%sM to %sM\nTm=%s"%(os.path.basename(datafile), str(round(minsteps/1000000.0,3)),str(round(maxtime/1000000.0,3)), str(round(Tmelt,4))))
		ax1.set_ylabel("Cv(T)")
		if Tunit=="Kelvin":
			ax1.set_xlabel("T(Kelvin)")
		else:
			ax1.set_xlabel("T(model)")
		
		ax1.plot(tkel,cvs,marker="x",label="Cv",linewidth=2.0)
		
		print "Tm=%f (i=%i); Cv(Tm)=%f"%(Tmelt,Tmelti, maxcv)
		
		if polyfitDegree != None:
			fitcvs, fitTmelt, fitTmelti, fitmaxcv, fineDCv  = getCv(temps, etotPerTemp, minTmCut=minTmCut, maxTmCut=maxTmCut, polyfitDegree=polyfitDegree, Tunit=Tunit)
			finegrid = list(numpy.arange(min(tkel),max(tkel),(max(tkel)-min(tkel)) / 100.0) )
			ax1.plot(finegrid,fitcvs,marker="",label="fitted Cv (deg %i)"%polyfitDegree,linewidth=2.0)
			#ax2 = ax1.twinx()
			#ax2.plot(finegrid,fineDCv,marker="",label="fitted Cv (deg %i)"%polyfitDegree,linewidth=1.0)
			#ax2.set_label("2nd derivative")
			print "Polynomial fit to Cv curve, degree=",polyfitDegree
			print "fitted Tm, max_Cv:",fitTmelt,fitmaxcv
			ax1.set_title(ax1.get_title()+" (polyfit(%i) Tm=%s)"%(polyfitDegree, str(round(fitTmelt,4))))
		print "Tm forced to be between %f and %f"%(minTmCut, maxTmCut)
		
		if vantHoff:
			ax1.plot(tkel,baselineCorrectedCv,marker="+",label="Cv (for coop.)",linewidth=2.0)
			ax1.plot(tkel,left_linCv,marker=".",label="left baseline",linewidth=2.0)
			ax1.plot(tkel,right_linCv,marker=".",label="right baseline",linewidth=2.0)
			ax1.plot([tkel[left_end_index],tkel[left_end_index]],[0,maxcv],label="T0",marker="", color="r",linewidth=1,linestyle=":")
			ax1.plot([tkel[right_begin_index],tkel[right_begin_index]],[0,maxcv],label="T1",marker="", color="r",linewidth=1,linestyle=":")
		
			g = [gaussfuncScalar(c, gauss_popt[0], gauss_popt[1],m=Tmelt) for c in baselineCorrectedCv]
			print g, "g"
			ax1.plot(tkel,g,marker="",label=" Gauss of Cv (for coop.)",linewidth=2.0)
		
		if optTmethod == "curve":
			ax1.plot(Tnew,newC,marker="o",label="Cv",linewidth=2.0)
		
		if tgrid_adjustRange:
			ax1.plot(tkel,[-200 for i in xrange(len(tkel))],marker="o",label="Told")
			ax1.plot(Tnew,[-100 for i in xrange(len(Tnew))],marker="v",label="Tnew")
			ax1.set_ylim(-300)
		
		if False: #experimental feature: Diffusivity
			ax2 = ax1.twinx()
			ax2.plot([str(tkel[i]+ (tkel[i+1] - tkel[i])/2.0) for i in range(0,len(tkel)-1)], Diff,  marker="", color="k",linewidth=1,linestyle=":")
			#ax2.plot(tkel,fT,marker="", color="k",linewidth=1,linestyle=":")
			#ax2.plot(tkel,fT,label="Cv",marker="o",color="b")
			ax2.set_ylabel("ln D(T)")
			#ax2.set_ylim(-0.0005)
		
		ax1.legend(loc=0)
		p.savefig(os.path.basename(datafile).replace("_rt.txt","_Cv.png"), dpi=figlowdpi)
		if showPlots: p.show()
		
		#write Cv data
		cvout = open(os.path.basename(datafile).replace("_rt.txt","_Cv.csv"),"w")
		cvout.write(",".join(["T","Cv"])+"\n")
		for i in xrange(len(tkel)):
			cvout.write(",".join(map(str,[tkel[i],cvs[i]]))+"\n")
		cvout.close()
		
		if polyfitDegree != None:
			#write fitted Cv data
			cvout = open(os.path.basename(datafile).replace("_rt.txt","_CvFit.csv"),"w")
			cvout.write(",".join(["T","fitted_Cv"])+"\n")
			for i in xrange(len(finegrid)):
				cvout.write(",".join(map(str,[finegrid[i],fitcvs[i]]))+"\n")
			cvout.close()
		#sys.exit()
		
	if trackStateTransitions:
		states = []
		#states.append([0.7,1.0,0.0,0.4]) # native A
		#states.append([0.0,0.7,0.0,0.4]) # unfolded
		#states.append([0.0,0.7,0.4,1.0]) # native B
		states.append([0.7,1.0,0.0,0.4]) # native A
		states.append([0.0,0.7,0.0,0.7]) # unfolded
		states.append([0.0,0.7,0.7,1.0]) # native B
		labels = ['A','U','B']
		#states.append([0.75,1.0,0.0,0.22]) # A
		#states.append([0.6,0.75,0.0,0.22]) # trA
		#states.append([0.0,0.6,0.0,0.22]) # U
		#states.append([0.0,0.6,0.22,0.4]) # UB
		#states.append([0.0,0.6,0.0,0.4]) # U
		#states.append([0.0,0.6,0.4,0.45]) # trB1
		#states.append([0.0,0.6,0.45,0.6]) # IB
		#states.append([0.0,0.6,0.6,0.66]) # trB2
		#states.append([0.0,0.6,0.66,1.0]) # B
		#labels = ['A','trA','U','UB','trB1','IB','trB2','B']
		#labels = ['A','trA','U','trB1','IB','trB2','B']
		
		tr = [[] for r in xrange(len(raw_data_per_replica))]
		n_trans = 0
		transitions_all = [[0 for s2 in xrange(len(states))] for s1 in xrange(len(states))]
		transitions_rep = [[[0 for s2 in xrange(len(states))] for s1 in xrange(len(states))] for r in xrange(len(raw_data_per_replica))]
		print len(raw_data_per_replica), "trajectories/replicas found"
		
		for r in xrange(len(raw_data_per_replica)):
			print "R",r," - ",len(raw_data_per_replica[r])
			for i in xrange(len(raw_data_per_replica[r])):
				#print raw_data_per_replica[r][i]
				q1 = raw_data_per_replica[r][i][3]
				q2 = raw_data_per_replica[r][i][4]
				for s in xrange(len(states)):
					if isInState(q1,q2,states[s]):
						
						if len(tr[r])==0:
							tr[r].append(s)
						else:
						#	if tr[r][-1] != s:
							transitions_rep[r][tr[r][-1]][s] += 1
							transitions_all[tr[r][-1]][s] += 1
							n_trans += 1
							tr[r].append(s)
								
			#print numpy.array(transitions_rep[r])
		print "\t"+"\t".join(labels)
		transitions_all = numpy.array(transitions_all)/float(n_trans)
		
		print "\n".join([labels[t1]+"\t"+"\t".join([str(round(transitions_all[t1][t2],4)) for t2 in xrange(len(transitions_all[t1]))]) for t1 in xrange(len(transitions_all))])
		transitions_all2 = numpy.linalg.matrix_power(numpy.array(transitions_all)/float(n_trans) , 2)
		print "A->B\t"+"\t".join([str(transitions_all[i-1][i]) for i in range(1,len(states))])
		print "A->A\t"+"\t".join([str(transitions_all[i][i]) for i in range(0,len(states))])
		transFile = open(os.path.basename(datafile).replace("_rt.txt",".csv"),'w')
		transFile.write(os.path.basename(datafile).replace("_rt.txt","")+"\t"+"\t".join(labels)+"\n")
		transFile.write(os.path.basename(datafile).replace("_rt.txt","_")+"A->B\t"+"\t".join([str(transitions_all[i-1][i]) for i in range(1,len(states))])+"\n")
		transFile.write(os.path.basename(datafile).replace("_rt.txt","_")+"A->A\t"+"\t".join([str(transitions_all[i][i]) for i in range(0,len(states))])+"\n")
		transFile.close()
		#print
		#print "\n".join(["\t".join([str(t2) for t2 in t1]) for t1 in transitions_all2])
		if True:
			p.matshow(transitions_all)
			for i in xrange(len(labels)):
				p.text(i,i,labels[i])
			#p.xlabels(labels)
			#p.ylabels(labels)
			p.colorbar()
			p.savefig(os.path.basename(datafile).replace("_rt.txt",".png"), dpi=figlowdpi)
		#sys.exit()
	if constTempStats:
		doConstTempStats(raw_data_per_replica)
	
	# HISTOGRAMS
	if not(print_free_energy or print_histogram or freeEnergy1D or freeEnergy2D): sys.exit()
	
	histokeys, histokeysA, histokeysB, countPerTemp, histoPerTemp3Deff, totalCount, histoPerTempQAE, histoPerTempQBE, histoPerTempEtot,histoPerTempQ1,histoPerTempQ2 = getHistograms(raw_data, QAbins, QBbins, Etotbins, temps)
	#print histoPerTemp2D
	
	# PRINT FREE ENERGY
	
	if print_free_energy or print_histogram:
		#print "Free energy surfaces:"
		sq = int(math.ceil(math.sqrt(len(temps))))
		fig = p.Figure(figsize=(12, 12),dpi=200)
		#print datafile
		p.figtext(0.5, 0.98, datafile, ha='center', color='black', weight='bold', size='large')
		p.figtext(0.5, 0.95, "time (min,max): %i,%i"%(minsteps,maxtime),ha='center', color='black', weight='bold', size='large')
	
		for t in xrange(len(temps)):
			ax = p.subplot(sq, sq, t+1)
			if Tunit=="Kelvin":
				ax.set_title("i = %i; beta = %s (%sK)"%(t, str(round(temps[t][2],2)), str(round(temps[t][1],2) ) ) )
			else:
				ax.set_title("i = %i; Tmodel = %s"%(t, str(str(round(temps[t][0],4) ) ) ))
			if t in [i for i in range(len(temps)-sq, len(temps) ) ]:	
				p.xlabel(QBlabel)
			else: 
				p.xlabel("")
		
			if t in [i for i in xrange(len(temps)) if i==0 or i%sq==0 ]: 
				p.ylabel(QAlabel)
			else: 
				p.ylabel("")
			#p.ylim(Pmin,Pmax)
			#p.colorbar(label="-T*ln(P)")#"-T*ln(P)"
			
			if print_free_energy:
			
				maximum = max([ max([ -temps[t][0] * math.log( j/float(countPerTemp[t]) ) if j != 0 else None for j in i]) for i in histoPerTemp2D[t] ])
			
				dat = []
				datB = []
				quadcount = [0,0,0,0]
			
				for i in xrange(len(histoPerTemp2D[t])):
					dat.append([])
					datB.append([])
					for j in xrange(len(histoPerTemp2D[t][i])):
				
						if histoPerTemp2D[t][i][j] != 0.0:
					
							dat[i].append(-temps[t][0] * math.log( histoPerTemp2D[t][i][j]/float(countPerTemp[t]) ) )
							datB[i].append(-temps[t][0] * math.log( histoPerTemp2D[t][i][j]/float(countPerTemp[t]) ) )
						else:
							dat[i].append(maximum)
							datB[i].append(0.0)
				#print countPerQuadrantPerT(t,[datB], QAbins, QBbins, qa_bisect, qb_bisect,verbose=False)[0], "T%i ="%t,temps[t][1]
			
				p.contourf(QBbins, QAbins, dat, ncontours)
				cbar = p.colorbar()
			
				if t in [i for i in xrange(len(temps)) if (i+1)%sq==0 ]: 
					cbar.set_label("-T*ln(P)")
			elif print_histogram:
				p.contourf(QBbins,QAbins, [ [ j/float(countPerTemp[t]) for j in i] for i in histoPerTemp2D[t] ] , ncontours)
		p.show()
	
	if finalTemps == []:
		if len(temps) > 1:
			cvs, Tmelt, Tmelti, maxcv, dCv  = getCv(temps, etotPerTemp,maxTmCut=maxTmCut,minTmCut=minTmCut,Tunit=Tunit)
			#acmeans = autoCorrelation(temps,etotPerTemp,plots=False)
			#minac = min(acmeans)
			#minaci = acmeans.index(minac)
			#Tnew, cvpeaks = optTbyOscillator(tkel, cvs, rt=5.0)
			#finalTemps = [tkel[i] for i in cvpeaks if not tkel[i] in finalTemps]
			finalTemps.append(tkel[Tmelti])
			#finalTemps.append(tkel[minaci])
			#if not freeEnergy1D:
			#	finalTemps.append(tkel[0])
			finalTemps = list(reversed(sorted(finalTemps)))
		else:
			finalTemps.append(tkel[0])
	
	print finalTemps
	
	if markStructuresFromFile:
		cfile = open(centroids_fname)
		abspaths = [os.readlink(os.path.abspath(i)) if os.path.islink(os.path.abspath(i)) else os.path.abspath(i) for i in cfile.readline().strip().split(",")]
		centroids = [os.path.basename(i) for i in abspaths]
		from kmeans import getAdjacencyMatrixFromFeature
		A = getAdjacencyMatrixFromFeature(abspaths, centroidDistCutoff, feature, ref_pdb=superpos_refpdbs, atom_types=superpos_atypes)
		#A = numpy.power(A,2)
		if markStructuresFromFile:
			centroidString = "%scentroids_%s_"%(len(centroids),feature)
		else:
			centroidString = ""
		if nativenessQ in ["True","true"]:
			confilename = os.path.basename(datafile).replace("_rt.txt","_%s_F_%s_%s_connectivities.dat"%(centroidString,QAlabel,QBlabel))
		else:
			confilename = os.path.basename(datafile).replace("_rt.txt","_%s_F_EA_EB_connectivities.dat"%(centroidString))
		Afile = open(confilename,"w")
		
		clusters = []
		for line in cfile.readlines():
			clusters.append(line.strip().split(","))
		assert len(clusters)==len(centroids)
		nMembersTotal = 0
		for c in clusters:
			nMembersTotal += len(c)
		minClusterSize = minClusterSizeFraction*nMembersTotal
		print "min cluster size:",minClusterSize
		#centroids = [os.path.basename(i.strip()) for i in open(centroids_fname).readlines()]
		pdb2data = {}
		p2d_file = open(pdb2data_fname)
		p2d_header = p2d_file.readline().strip().split(",")
		for line in p2d_file.readlines():
			tmp = line.strip().split(",")
			pdb2data[os.path.basename(tmp[0])] = [ float(i) for i in tmp[1:] ]
		for c in centroids:
			assert c in pdb2data,c
		
		clusterNeighbours = {}
		adjPairs = []
		connectivities={}
		for i in xrange(len(centroids)):
			connectivities[os.path.basename(centroids[i])] = 0
			Afile.write( ",".join(map(str,A[i]))+"\n" )
			clusterNeighbours[i] = {}
			for j in xrange(len(centroids)):
				
				if A[i][j] != 0.0:
					connectivities[os.path.basename(centroids[i])] += 1
					if A[i][j] < centroidDistCutoff and i != j:
						clusterNeighbours[i][j] = A[i][j]
					if i>j:
						
						adjPairs.append( (os.path.basename(centroids[i]),os.path.basename(centroids[j]),A[i][j] ))
			#print i,connectivities[os.path.basename(centroids[i])]
	
		Afile.close()
	
	# WHAM
	if freeEnergy1D:
		# histo, energy, q
		if True:
			for t in xrange(len(temps)):
				print t,temps[t]
			if constTempStats:
				Tmelti = 0
			
			if Tmelti != None:
				raw = raw_input("Choose temperature index [Enter:%i]:"%Tmelti)
				if raw.strip() != "":
					Tmelti = int(raw)
			else:
				Tmelti = int(raw_input("Choose temperature index:"))
			p.clf()
			histo_energy = [[ histoPerTempEtot[t][e]  for e in xrange(nEtotbins)] for t in xrange(len(temps))]
			p.plot(Etotbins, histo_energy[Tmelti])
			p.xlabel("energy")
			p.ylabel("count")
			p.savefig(os.path.basename(datafile).replace("_rt.txt","_histE.png"), dpi=200)
			if showPlots: p.show()
			
			p.clf()
			histo_q1 = [[ histoPerTempQ1[t][q]  for q in xrange(nQAbins)] for t in xrange(len(temps))]
			p.plot(QAbins, histo_q1[Tmelti])
			p.xlabel(QAlabel)
			p.ylabel("count")
			p.savefig(os.path.basename(datafile).replace("_rt.txt","_histQA.png"), dpi=200)
			if showPlots: p.show()
			
			p.clf()
			histo_q2 = [[ histoPerTempQ2[t][q]  for q in xrange(nQBbins)] for t in xrange(len(temps))]
			p.plot(QBbins, histo_q2[Tmelti])
			p.xlabel(QBlabel)
			p.ylabel("count")
			p.savefig(os.path.basename(datafile).replace("_rt.txt","_histQB.png"), dpi=200)
			if showPlots: p.show()
	
			p.clf()
			histo_E_q1 = [[ [ histoPerTempQAE[t][(q,e)] if (q,e) in histoPerTempQAE[t] else 0.0 for e in xrange(nEtotbins)]   for q in xrange(nQAbins)] for t in xrange(len(temps))]
			p.contourf( Etotbins,QAbins, histo_E_q1[Tmelti] , ncontours)
			p.colorbar()
			p.xlabel("E")
			p.ylabel(QAlabel)
			p.savefig(os.path.basename(datafile).replace("_rt.txt","_histEQA.png"), dpi=200)
			if showPlots: p.show()
			
			p.clf()
			histo_E_q2 = [[ [ histoPerTempQBE[t][(q,e)] if (q,e) in histoPerTempQBE[t] else 0.0 for e in xrange(nEtotbins)]   for q in xrange(nQBbins)] for t in xrange(len(temps))]
			p.contourf( Etotbins,QBbins, histo_E_q2[Tmelti] , ncontours)
			p.colorbar()
			p.xlabel("E")
			p.ylabel(QBlabel)
			p.savefig(os.path.basename(datafile).replace("_rt.txt","_histEQB.png"), dpi=200)
			if showPlots: p.show()
		
		#QA
		if True:
			P_unb_norm = wham(finalTemps, histoPerTempQAE, countPerTemp, Emidpoints, max_wham_iterations, histokeysA, temps, useDecimalPackage=useDecimalPackage, verbose=False, wham_cut=1e-3, Tunit="Profasi" )
			#wham QA 1D profile
			if True:
				sq = int(math.ceil(math.sqrt(len(finalTemps))))
				p.clf()
				fig = p.Figure(figsize=(12, 12),dpi=200)
				title = os.path.basename(datafile).replace("_rt.txt","")+";T%i-%i;%ixWHAM;P>%s"%(minsteps,maxtime,max_wham_iterations,str(whamMinP))
				p.figtext(0.5, 0.965, title,ha='center', color='black', weight='bold', size='large')

				for t in xrange(len(finalTemps)):
					if P_unb_norm[t] == None: continue
					ax = p.subplot(sq, sq, t+1)
					ax.set_title("i = %i; T = %s "%(t, str(round(finalTemps[t],4) ) ) )
					if t in [i for i in range(len(finalTemps)-sq, len(finalTemps) ) ]:
						p.xlabel(QAlabel)
					else: 
						p.xlabel("")

					if t in [i for i in xrange(len(finalTemps)) if i==0 or i%sq==0 ]: 
						p.ylabel("free Energy")
					else: 
						p.ylabel("")
					#p.ylim(Pmin,Pmax)

					summedEbins = [ sum([ P_unb_norm[t][(q,e)] if (q,e) in P_unb_norm[t] else 0.0 for e in xrange(nEtotbins)] )  for q in xrange(nQAbins)]
					#print summedEbins
					
					if Tunit=="Kelvin":
						prfT = 1.0/KtoPRFbeta(finalTemps[t])
					else:
						prfT = finalTemps[t]
					
					minimum =  min([ -prfT * math.log(summedEbins[q]) if summedEbins[q] > 1e-50 else 1e100  for q in xrange(nQAbins)])
					
					dat =  [ -prfT * math.log(summedEbins[q])-minimum if summedEbins[q] > 1e-50 else None  for q in xrange(nQAbins)]
					
					print "QA(%s),F(QA;%s;Tm=%i)"%(datafile,datafile,finalTemps[t])
					for i in xrange(len(dat)):
						print ",".join([str(QAmidpoints[i]),str(dat[i])])
					
					p.plot(QAbins,dat)
					
					if False:
						minFold, minbin_Fold = findMinimumInHisto1D(dat, QAbins, 0.6, 1.0)
						minUnfold, minbin_Unfold = findMinimumInHisto1D(dat, QAbins, 0.0, 0.6)
						
						print "Temperature",finalTemps[t]
						print "Free Energy Minimum at Folded for Q bin (>0.6):",minFold, QAbins[minbin_Fold]
						print "Free Energy Minimum at Unfolded for Q bin (<0.6):",minUnfold, QAbins[minbin_Unfold]
				p.savefig(os.path.basename(datafile).replace("_rt.txt","_QA1D.png"))
				if showPlots: p.show()
				p.clf()
		
			#free energy E,QA  2D
			if True:
				sq = int(math.ceil(math.sqrt(len(finalTemps))))
				p.clf()
				
				fig = p.Figure(figsize=(12, 12),dpi=200)
				title = os.path.basename(datafile).replace("_rt.txt","")+";T%i-%i;%ixWHAM;P>%s"%(minsteps,maxtime,max_wham_iterations,str(whamMinP))
				p.figtext(0.5, 0.965, title,ha='center', color='black', weight='bold', size='large')

				for t in xrange(len(finalTemps)):
					if P_unb_norm[t] == None: continue
					ax = p.subplot(sq, sq, t+1)
					ax.set_title("i = %i; T = %s "%(t, str(round(finalTemps[t],4) ) ) )
					if t in [i for i in range(len(finalTemps)-sq, len(finalTemps) ) ]:
						p.xlabel(QAlabel)
					else: 
						p.xlabel("")

					if t in [i for i in xrange(len(finalTemps)) if i==0 or i%sq==0 ]: 
						p.ylabel("Total Energy")
					else: 
						p.ylabel("")
					
					if Tunit=="Kelvin":
						prfT = 1.0/KtoPRFbeta(finalTemps[t])
					else:
						prfT = finalTemps[t]
					
					try:
						minimum =  min([min([ -prfT * math.log(P_unb_norm[t][(q,e)])  for q in xrange(nQAbins) if (q,e) in P_unb_norm[t]]) for e in xrange(nEtotbins)])
			
						maximum =  max([max([ -prfT * math.log(P_unb_norm[t][(q,e)])  for q in xrange(nQAbins) if (q,e) in P_unb_norm[t]]) for e in xrange(nEtotbins)])
					except ValueError as ve:
						print ve
						minimum = 0.0
						maximum = 1.0
					dat =  [[ -prfT * math.log(P_unb_norm[t][(q,e)])-minimum if (q,e) in P_unb_norm[t] else  maximum-minimum  for q in xrange(nQAbins)] for e in xrange(nEtotbins)]
					#print dat
					if False:
						minQA, minbinA_QA, minbinB_QA = findMinimumInHisto(dat, Etotbins, QAbins, -200, 100, 0.8, 1.0)
						minQB, minbinA_QB, minbinB_QB = findMinimumInHisto(dat, Etotbins, QAbins, -200, 100, 0.2, 0.4)
			
						print "free energy minimum between QA=0.8 and 1.0: ",minQA, Etotbins[minbinA_QA], QAbins[minbinB_QA]
						print "free energy minimum between QA=0.2 and 0.4:",minQB, Etotbins[minbinA_QB], QAbins[minbinB_QB]
			
					#dat =  [histoPerTempQA[t][q]   for q in xrange(nQAbins)]
					levels = [0.0+i*(8.0/ncontours) for i in xrange(ncontours)]
					p.contourf(QAbins, Etotbins, dat , ncontours,levels=levels,extend="max")
					cbar = p.colorbar()
			
				p.title("Free Energy")
				p.savefig(os.path.basename(datafile).replace("_rt.txt","_EQA2D.png"))
				if showPlots: p.show()
		
		#QB
		if True:
			P_unb_norm = wham(finalTemps, histoPerTempQBE, countPerTemp, Emidpoints, max_wham_iterations, histokeysB, temps, useDecimalPackage=False, verbose=False, wham_cut=1e-3, Tunit="Profasi")
			#wham QB 1D profile
			if True:
				sq = int(math.ceil(math.sqrt(len(finalTemps))))
				p.clf()
				
				fig = p.Figure(figsize=(12, 12),dpi=200)
				title = os.path.basename(datafile).replace("_rt.txt","")+";T%i-%i;%ixWHAM;P>%s"%(minsteps,maxtime,max_wham_iterations,str(whamMinP))
				p.figtext(0.5, 0.965, title,ha='center', color='black', weight='bold', size='large')	
				for t in xrange(len(finalTemps)):
					if P_unb_norm[t] == None: continue
					ax = p.subplot(sq, sq, t+1)
					ax.set_title("i = %i; T = %s "%(t, str(round(finalTemps[t],4) ) ) )
					if t in [i for i in range(len(finalTemps)-sq, len(finalTemps) ) ]:
						p.xlabel(QBlabel)
					else: 
						p.xlabel("")

					if t in [i for i in xrange(len(finalTemps)) if i==0 or i%sq==0 ]: 
						p.ylabel("free Energy")
					else: 
						p.ylabel("")
					#p.ylim(Pmin,Pmax)

					summedEbins = [ sum([ P_unb_norm[t][(q,e)] if (q,e) in P_unb_norm[t] else 0.0 for e in xrange(nEtotbins)] )  for q in xrange(nQBbins)]
					#print summedEbins
					if Tunit=="Kelvin":
						prfT = 1.0/KtoPRFbeta(finalTemps[t])
					else:
						prfT = finalTemps[t]
						
					minimum =  min([ -prfT * math.log(summedEbins[q]) if summedEbins[q] > 1e-50 else 1e100  for q in xrange(nQBbins)])

					dat =  [ -prfT * math.log(summedEbins[q])-minimum if summedEbins[q] > 1e-50 else None  for q in xrange(nQBbins)]
					#print dat
					print "QB(%s),F(QB;%s;Tm=%f)"%(datafile,datafile,finalTemps[t])
					for i in xrange(len(dat)):
						print ",".join([str(QBmidpoints[i]),str(dat[i])])
		
					p.plot(QBbins,dat)
			
				p.savefig(os.path.basename(datafile).replace("_rt.txt","_QB1D.png"))
				if showPlots: p.show()
		
			#free energy E,QB  2D
			if True:
				sq = int(math.ceil(math.sqrt(len(finalTemps))))
				p.clf()
				
				fig = p.Figure(figsize=(12, 12),dpi=200)
				title = os.path.basename(datafile).replace("_rt.txt","")+";T%i-%i;%ixWHAM;P>%s"%(minsteps,maxtime,max_wham_iterations,str(whamMinP))
				p.figtext(0.5, 0.965, title,ha='center', color='black', weight='bold', size='large')

				for t in xrange(len(finalTemps)):
					if P_unb_norm[t] == None: continue
					ax = p.subplot(sq, sq, t+1)
					ax.set_title("i = %i; T = %s "%(t, str(round(finalTemps[t],4) ) ) )
					if t in [i for i in range(len(finalTemps)-sq, len(finalTemps) ) ]:
						p.xlabel(QBlabel)
					else: 
						p.xlabel("")

					if t in [i for i in xrange(len(finalTemps)) if i==0 or i%sq==0 ]: 
						p.ylabel("Total Energy")
					else: 
						p.ylabel("")
					
					if Tunit=="Kelvin":
						prfT = 1.0/KtoPRFbeta(finalTemps[t])
					else:
						prfT = finalTemps[t]
					#summedEbins = [ sum([ P_unb_norm[t][(q,e)] if (q,e) in P_unb_norm[t] else 0.0 for e in xrange(nEtotbins)] )  for q in xrange(nQAbins)]
					#print summedEbins
					minimum =  min([min([ -prfT * math.log(P_unb_norm[t][(q,e)])  for q in xrange(nQBbins) if (q,e) in P_unb_norm[t]]) for e in xrange(nEtotbins)])
		
					maximum =  max([max([ -prfT * math.log(P_unb_norm[t][(q,e)])  for q in xrange(nQBbins) if (q,e) in P_unb_norm[t]]) for e in xrange(nEtotbins)])
		
					dat =  [[ -prfT * math.log(P_unb_norm[t][(q,e)])-minimum if (q,e) in P_unb_norm[t] else  maximum-minimum  for q in xrange(nQBbins)] for e in xrange(nEtotbins)]
					#print dat
	
					#dat =  [histoPerTempQA[t][q]   for q in xrange(nQAbins)]
		
					if False:
						minQA, minbinA_QA, minbinB_QA = findMinimumInHisto(dat, Etotbins, QBbins, -200, 100, 0.7, 1.0)
						minQB, minbinA_QB, minbinB_QB = findMinimumInHisto(dat, Etotbins, QBbins, -200, 100, 0.1, 0.3)
						print "free energy minimum between QB=0.7 and 1.0: ",minQA, Etotbins[minbinA_QA], QBbins[minbinB_QA]
						print "free energy minimum between QB=0.1 and 0.3: ",minQB, Etotbins[minbinA_QB], QBbins[minbinB_QB]
		
					p.contourf(QBbins, Etotbins, dat , ncontours)
					cbar = p.colorbar()
		
				p.title("Free Energy")
				p.savefig(os.path.basename(datafile).replace("_rt.txt","_EQB2D.png"))
				if showPlots: p.show()
				p.clf()
	if freeEnergy2D:
		p.clf()
		p.close()
		P_unb_norm = wham(finalTemps, histoPerTemp3Deff, countPerTemp, Emidpoints, max_wham_iterations, histokeys, temps, useDecimalPackage=False, verbose=False, wham_cut=1e-3, Tunit="Profasi")
		# plot  2D surface wham
		if True:
			sq = int(math.ceil(math.sqrt(len(finalTemps))))
			fig = p.Figure(figsize=(15, 15),dpi=800)
			bg = fig.patch
			bg.set_facecolor(BGcolor)
			title = os.path.basename(datafile).replace("_rt.txt","")+";T%i-%i"%(minsteps,maxtime)
			
			if drawMainTitle:
				p.figtext(0.5, 0.965, title,ha='center', color='black', weight='normal', size='medium')

			for t in xrange(len(finalTemps)):
				if P_unb_norm[t] == None: continue
				ax = p.subplot(sq, sq, t+1)
				ax.set_axis_bgcolor(BGcolor)
				ax.set_aspect(aspect)
				if not showAxes:
					ax.set_axis_off()
				
				if drawSubTitle:
					if Tunit=="Kelvin":
						ax.set_title("i = %i; beta = %s (%sK)"%(t, str(round(KtoPRFbeta(finalTemps[t]),2)), str(round(finalTemps[t],2) ) ) )
					else:
						ax.set_title("i = %i; Tmodel = %s "%(t, str(round(finalTemps[t],4) ) ) )
				if t in [i for i in range(len(finalTemps)-sq, len(finalTemps) ) ]:	
					p.xlabel(QBlabel)
				else: 
					p.xlabel("")

				if t in [i for i in xrange(len(finalTemps)) if i==0 or i%sq==0 ]: 
					p.ylabel(QAlabel)
				else: 
					p.ylabel("")
				#p.ylim(Pmin,Pmax)
				
				#print P_unb_norm[t]
				summedEbins = [ [ sum([ P_unb_norm[t][(q1,q2,e)] if (q1,q2,e) in P_unb_norm[t] else 0.0 for e in xrange(nEtotbins)] )  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
				#print summedEbins
				
				if Tunit=="Kelvin":
					prfT = 1.0/KtoPRFbeta(finalTemps[t])
				else:
					prfT = finalTemps[t]
				
				minimum = min([ min([ -prfT * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else 1e900 for q2 in xrange(nQBbins) ]) for q1 in xrange(nQAbins)])
				maximum = max([ max([ -prfT * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else 1e-900 for q2 in xrange(nQBbins) ]) for q1 in xrange(nQAbins)])
				
				print minimum,maximum
				if maxFreeEnergy == None:
					dat = [ [ -prfT * math.log(summedEbins[q1][q2])-minimum if summedEbins[q1][q2] != 0.0 else maximum-minimum  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
				else:
					dat = [ [ -prfT * math.log(summedEbins[q1][q2])-minimum if summedEbins[q1][q2] != 0.0 else maxFreeEnergy  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
				#print dat
				
				#datB = [ [ -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else 0.0  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
				
				if exportFmatrix:
					if nativenessQ in ["True","true"]:
						Fname = os.path.basename(datafile).replace("_rt.txt","_F_%s_%s.dat"%(QAlabel,QBlabel))
					else:
						Fname = os.path.basename(datafile).replace("_rt.txt","_F_EA_EB.dat")
					print Fname
					outF = open(Fname,"w")
					
					nanDat = dat
					#nanDat = [ [ -(1.0/KtoPRFbeta(finalTemps[t])) * math.log(summedEbins[q1][q2]) if summedEbins[q1][q2] != 0.0 else None  for q2 in xrange(nQBbins)] for q1 in xrange(nQAbins)]
					
					for i in nanDat:
						outF.write(",".join(map(str,i))+"\n")
					outF.close()
				
				minQA, minbinA_QA, minbinB_QA = findMinimumInHisto(dat, QAbins, QBbins, qa_bisect, max(1.0,EAmin), 0.0, qb_bisect)
				minQB, minbinA_QB, minbinB_QB = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, qb_bisect, max(1.0,EBmin))
				minUnf, minbinA_Unf, minbinB_Unf = findMinimumInHisto(dat, QAbins, QBbins, 0.0, qa_bisect, 0.0, qb_bisect)
				
				print "Temperature",finalTemps[t]
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
					if markStructuresFromFile:
						
						p.contour(QBbins, QAbins, dat, ncontours, extend="max", alpha=0.5,linewidths=1.0,cmap="Dark2",figure=fig, zorder=2)
					else:
						p.contourf(QBbins, QAbins, dat, ncontours,levels=levels, extend="max",figure=fig, zorder=2)
				else:
					if markStructuresFromFile:
						p.contour(QBbins, QAbins, dat , ncontours, alpha=0.75,linewidths=0.5,cmap="Dark2",figure=fig, zorder=2)
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
				
					if t in [i for i in xrange(len(temps)) if (i+1)%sq==0 ]: 
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
			
			if markStructuresFromFile:
				
				cluster_histo_in_range = [{} for i in xrange(len(drawRectangles))]
				markercolors = []
				if analyseRectangles and exportPDBsInRange:
					for r in drawRectangles:
						if markStructuresFromFile:
							centroidString = "%scentroids_%s_"%(len(centroids),feature)
						else:
							centroidString = ""
						if nativenessQ in ["True","true"]:
							foldername = os.path.basename(datafile).replace("_rt.txt","_%s_F_%s_%s_%s_%s_%s_%s"%(centroidString, QAlabel, QBlabel , str(r[0]), str(r[1]), str(r[2]), str(r[3])))
						else:
							foldername = os.path.basename(datafile).replace("_rt.txt","_%s_F_EA_EB_%s_%s_%s_%s"%(centroidString , str(r[0]), str(r[1]), str(r[2]), str(r[3])))
						if not os.path.exists(foldername):
							os.mkdir(foldername)
				
				if drawRectangles != []:
					#p = []
					for rect in drawRectangles:
						print rect
						ax.add_patch(mpatches.Rectangle((rect[0],rect[2]),rect[1]-rect[0],rect[3]-rect[2], facecolor=rect[4], fill=True, alpha=rect[4][3], edgecolor="none" ))
						#p.append(  )
					#collection = PatchCollection(patches)#, cmap=cm.hsv, alpha=0.3)
					#collection.set_array(numpy.array(colors))
					#ax.add_collection(collection)
					
				
				for ci in xrange(len(centroids)):
					qas=[]
					qbs=[]
					#print ci, len(clusters[ci])
					for i in xrange(len(clusters[ci])):
						if os.path.islink(clusters[ci][i]):
							cl = os.path.basename(os.readlink(clusters[ci][i]))
						else:
							cl = os.path.basename(os.path.abspath(clusters[ci][i]))
						tmp = pdb2data[cl][1:]
						
						#print tmp
						if nativenessQ in ["True","true","1"]:
							QA = float(tmp[indexA])
							QB = float(tmp[indexB])
							EA=None
							EB=None
						else:	
							EA = float(tmp[indexA])
							EB = float(tmp[indexB])
							QA = -float(EA*weightA)/(EAmin*weightA)
							QB = -float(EB)/EBmin
							if QA < 0.0: QA = 0.0
							if QB < 0.0: QB = 0.0
						QA_bin = getBinIndex(QAbins, QA)
						QB_bin = getBinIndex(QBbins, QB)
						QA_mid = QAmidpoints[QA_bin]
						QB_mid = QBmidpoints[QB_bin]
						#print indexA,indexB,QA,QB,EA,EB,QA_bin,QB_bin,QA_mid,QB_mid
						#print ci+minClusterIndex,c,QA,QB
						
						if analyseRectangles: #QA > qa_cluster_range[0] and QA < qa_cluster_range[1] and QB > qb_cluster_range[0] and QB < qb_cluster_range[1]):
							inR = isPointInRectangles(QB,QA,drawRectangles)
							if inR >= 0:
								if onlyDrawMembersInRectangles:
									qas.append(QA)
									qbs.append(QB)
							
								if exportPDBsInRange:
									shutil.copy(clusters[ci][i], os.path.join(foldername, cl))
								
								if ci in cluster_histo_in_range[inR]:
									cluster_histo_in_range[inR][ci] += 1
								else:
									cluster_histo_in_range[inR][ci] = 1
							else:
								if not onlyDrawMembersInRectangles:
									qas.append(QA)
									qbs.append(QB)
						else:
							qas.append(QA)
							qbs.append(QB)	
									
					r = lambda: random.uniform(0,1)
					if drawCentroidRandomColors:
						col = (r(),r(),r(),1.0)
					else:
						col = uniformCentroidColor
					
					if drawClusterMembers:
						
						p.scatter(qbs,qas,s=0.1,c=col,edgecolor=col,facecolor=col,marker="+",figure=fig, zorder=3)
					
					markercolors.append(col)
				
				#print cluster_histo_in_range
				inRange_sum = 0
				for r in xrange(len(drawRectangles)):
					print drawRectangles[r][5]
					sortedClustersInRect = sorted(cluster_histo_in_range[r].items(), key=operator.itemgetter(1), reverse=True)
					for jj in xrange(highlightNtopClusters):
						highlightclusters.append(sortedClustersInRect[jj][0])
					for kk in sortedClustersInRect:
						print "c%i: %i"%(kk[0],kk[1])
						inRange_sum += kk[1]
					print inRange_sum, "conformations in box"
				
				if drawCentroidConnections:
					for a in adjPairs:
						p1 = a[0]
						p2 = a[1]
						if connectivities[p1] < minCentroidConn or connectivities[p2] < minCentroidConn:
							continue
						v = a[2]
						qa1,qb1 = getCoordPDB(p1,pdb2data,indexA,indexB)
						qa2,qb2 = getCoordPDB(p2,pdb2data,indexA,indexB)
						p.plot([qb1,qb2], [qa1,qa2], color='black',marker="",linewidth=v*connectStrength,linestyle="-",figure=fig, zorder=4)
				
				print "Cluster centroids:"
				for ci in xrange(len(centroids)):
					c = centroids[ci]
					
					tmp = pdb2data[c][1:]
					#print tmp
					if nativenessQ in ["True","true","1"]:
						QA = float(tmp[indexA])
						QB = float(tmp[indexB])
						EA=None
						EB=None
					else:	
						EA = float(tmp[indexA])
						EB = float(tmp[indexB])
						QA = -float(EA*weightA)/(EAmin*weightA)
						QB = -float(EB)/EBmin
						if QA < 0.0: QA = 0.0
						if QB < 0.0: QB = 0.0
					QA_bin = getBinIndex(QAbins, QA)
					QB_bin = getBinIndex(QBbins, QB)
					QA_mid = QAmidpoints[QA_bin]
					QB_mid = QBmidpoints[QB_bin]
					#print indexA,indexB,QA,QB,EA,EB,QA_bin,QB_bin,QA_mid,QB_mid
					
					nmembers = len(clusters[ci])
					size = 4+math.log(nmembers)
					
					
					print ci+minClusterIndex,"\t",c,"\t",round(QA,2),"\t",round(QB,2),"\t",nmembers,"\t",connectivities[os.path.basename(c)], "\t",[s[0] for s in sorted(clusterNeighbours[ci].items(), key=operator.itemgetter(1), reverse=True) ]
					#p.plot([QB],[QA],marker="x",color="1.0",markersize=10, markeredgewidth=3)
					#p.text(QB, QA, str(ci+minClusterIndex), fontsize=12)
					cen_QA = QA
					cen_QB = QB
					col=markercolors[ci]
					#print col,numpy.average(col)
					
					if drawCentroidCircles:
						#if connectivities[os.path.basename(c)] >= minCentroidConn and nmembers>=minClusterSize:
						if nmembers>=minClusterSize:
							if ci in highlightclusters:
								ecol = cluster_highlight_color
							else:
								ecol = "white"
							p.scatter([cen_QB],[cen_QA],marker="o",c=col,s=sizeboost*nmembers*0.1, edgecolor=ecol,label=str(ci+minClusterIndex),alpha=1.0,figure=fig, zorder=5)
							
					if drawCentroidLabels:
						#if connectivities[os.path.basename(c)] >= minCentroidConn:
						#if nmembers>=minClusterSize:
						if uniformCentroidColor:
							textcol = (0.0,0.0,0.0,1.0)
						else:
							textcol = col
						text = p.text(cen_QB+0.01, cen_QA+0.01, str(ci+minClusterIndex), fontsize=sizeboost*size*0.5, color=textcol, backgroundcolor=(1,1,1,0.25),figure=fig, zorder=6)
					#text.set_bbox(dict(facecolor="red", alpha=0.5, size=1))
					#print text.get_window_extent()
					#print text.get_bbox()
					#print
						#p.text(QB, QA, str(ci+minClusterIndex), fontsize=12)
			if drawCentroidLegend:
				p.legend(loc=0,ncol=4,fontsize=4,markerscale=0.25)
			if markStructuresFromFile:
				centroidString = "%scentroids_%s_"%(len(centroids),feature)
			else:
				centroidString = ""
			
			
			if nativenessQ in ["True","true"]:
				figname = os.path.basename(datafile).replace("_rt.txt","_%s_F_%s_%s.png"%(centroidString,QAlabel,QBlabel))
			else:
				figname = os.path.basename(datafile).replace("_rt.txt","_%s_F_EA_EB.png"%(centroidString))
			print figname
			
			
			if not showAxes:
				p.axis('off')
			
			p.savefig(figname, dpi=figdpi,facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight', pad_inches = 0)	
			if showPlots: p.show()

if __name__=="__main__":
	
	##################
	# HANDLE OPTIONS (-h for help)
	description="%prog "
	version="ProFSBM v1.0; written by Tobias Sikosek in 2013-2016 at the University of Toronto ( tsikosek@icloud.com )"
	usage="\n\tpython %prog [OPTIONS] data.csv"
	
	op=optparse.OptionParser(description=description, usage=usage, version=version)
	
	op.add_argument('-i','--input', action="store", type="inputfilename", dest="inputfilename",help="A CSV data file created by 'exTraj.py' from PROFASI simulation folder.")
	
	run_anatraj(opt, args)
