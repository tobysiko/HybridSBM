import sys

from libprofsbm import doKmeans4Feature, filterSymLinksBySource, splitModels, joinModels

def run_kmeans(infolder, k , feature, filters, ref_pdb):
	assert feature in ["superdihedralcontacts","contacts","dihedrals","superposition","rmsd"]
	
	if feature in ["superdihedralcontacts","superposition","rmsd"]:
		atom_types = ["CA","CB"]
		
		if "," in ref_pdb:
			tmp_pdb = ref_pdb.split(",")
		else:
			tmp_pdb = [ref_pdb]
		
		ref_pdb = []
		print tmp_pdb
		for p in tmp_pdb:
			ref_pdb.extend(splitModels(p))
		
		print ref_pdb
		if feature == "rmsd":
			tmp_pdb = "_".join(ref_pdb)
			joinModels(ref_pdb, tmp_pdb)
			ref_pdb = tmp_pdb
	else:
		ref_pdb = None
		atom_types = None
	
	pdbclusters, centroids = doKmeans4Feature(infolder, k, feature,njobs=4, alwaysOverwrite=False, verbose=True, ref_pdb=ref_pdb, atom_types=atom_types)

	#subfolders, cent = createLinkSubfolders(infolder, centroids, feature, pdbclusters, shortlinks=True)

	totalconf = sum([len(pdbclusters[i]) for i in xrange(k)])
	print totalconf, "conformations in total"
	
	histo = [[0 for j in filters] for i in xrange(k)] 

	for i in xrange(k):
		#sf = subfolders[i]
		for f in xrange(len(filters)):
			histo[i][f] = len(filterSymLinksBySource(pdbclusters[i],filters[f]))
		#print i,histo[i]

	print filters
	total_counts = [0 for i in filters]
	for h in histo:
		for i in xrange(len(h)):
			total_counts[i] += h[i]
	#print total_counts,sum(total_counts)

	histofile = open(infolder+"kmeans%i_%s_histo.txt"%(k,feature),'w')
	histofile.write("labels,"+",".join(filters)+"\n")
	histofile.write("count,"+",".join(map(str,total_counts))+"\n")

	for h in histo:
		norm = [ h[i]/float(totalconf) for i in xrange(len(filters)) ]
		histofile.write("%s,"%str(histo.index(h))+",".join(map(str,norm))+"\n")
		print histo.index(h), h, norm, sum(norm)
	histofile.close()

if __name__=="__main__":
	infolder = sys.argv[1]
	k =    int(sys.argv[2])
	feature =  sys.argv[3]
	ref_pdb = sys.argv[4]
	
	filters= ["GA","GB"]
	
	run_kmeans(infolder, k, feature, filters, ref_pdb)


