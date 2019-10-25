import subprocess,sys,glob,os.path,shutil,random,math


##########################

folders = sys.argv[1:]

shortLinkNames = False

skipClustering = True

targetNoClusters = 12

targetTotalPdb = 40000

histoLabels = ["GA","GB"]

pdb2data = {}

if len(folders) > 1:
	pdbfolder = mergeFoldersLinks(folders,n=targetTotalPdb/len(folders))
	for f in folders:
		h = addFile2Data(pdb2data, "%s_pdb2data.dat"%f)
	writeData2Pdb(pdb2data, h, "%s_pdb2data.dat"%pdbfolder)
else:
	pdbfolder = folders[0]
	h = addFile2Data(pdb2data, "%s_pdb2data.dat"%pdbfolder)


print pdbfolder
print len(pdb2data), "entries in pdb2data"


listfile, nfiles = getList4Folder(pdbfolder, n=targetTotalPdb)

print nfiles," pdb files found"


maxClusterSize = nfiles / targetNoClusters

print maxClusterSize,"=max cluster size"
print targetNoClusters,"clusters attempted"

options = "-is %i"%maxClusterSize

if skipClustering:
	print "You just wanted to merge folders. No clustering today... abort."
	sys.exit()
maxcluster_text = maxcluster(pdbfolder, listfile, options=options, mcexepath="maxcluster")

print "wrote maxcluster output: ",maxcluster_text

#maxcluster_text = sys.argv[1]
centroids, clusterMembers, id2pdb, histo = parseMaxclusterOutput(maxcluster_text, labels=histoLabels)

print histoLabels
total_counts = [0 for i in histoLabels]
for h in sorted(histo):
	for i in xrange(len(histo[h])):
		total_counts[i] += histo[h][i]
print total_counts

histofile = open(pdbfolder+"".join(options.replace("-","_").split())+"_histo.txt",'w')
histofile.write("labels,"+",".join(histoLabels)+"\n")
histofile.write("count,"+",".join(map(str,total_counts))+"\n")

for h in sorted(histo):
	norm = [ histo[h][i]/float(total_counts[i]) for i in xrange(len(histo[h])) ]
	histofile.write("%s,"%str(h)+",".join(map(str,norm))+"\n")
	print h, norm
histofile.close()

createLinkSubfolders(pdbfolder, options, centroids, clusterMembers, id2pdb, shortlinks=shortLinkNames)



