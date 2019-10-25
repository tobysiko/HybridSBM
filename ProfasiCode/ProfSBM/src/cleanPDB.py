
import sys
from libprofsbm import parseATOMline, AA_3_to_1, writeATOMdict2Line

testPDB="""ATOM    193  N   GLN A  27      29.724  31.279  24.636  1.00 15.04           N  
ATOM    194  CA  GLN A  27      28.470  31.644  25.313  1.00 16.50           C  
ATOM    195  C   GLN A  27      28.439  32.914  26.073  1.00 18.25           C  
ATOM    196  O   GLN A  27      27.478  33.680  26.054  1.00 21.76           O  
ATOM    197  CB  GLN A  27      27.876  30.468  26.136  1.00 20.07           C  
ATOM    198  CG AGLN A  27      27.570  29.232  25.290  0.50 12.45           C  
ATOM    199  CG BGLN A  27      26.388  30.644  26.494  0.50 28.90           C  
ATOM    200  CD AGLN A  27      26.956  28.056  26.052  0.50 20.66           C  
ATOM    201  CD BGLN A  27      26.083  30.009  27.823  0.50 30.43           C  
ATOM    202  OE1AGLN A  27      27.168  26.882  25.699  0.50 22.61           O  
ATOM    203  OE1BGLN A  27      25.473  30.610  28.715  0.50 34.69           O  
ATOM    204  NE2AGLN A  27      26.206  28.365  27.097  0.50 23.95           N  
ATOM    205  NE2BGLN A  27      26.508  28.766  27.938  0.50 41.63           N  """

#AA_3_to_1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}



def run_cleanPDB(arguments):
	if not type(arguments)==list: arguments = [arguments]
	TEST=False
	
	
	if not TEST:
		name = arguments[0]
		#print name
		outpdb = open(name[:-4]+"_clean.pdb", 'w')
		lines = open(name).readlines()
	else:
		print testPDB
		print
		lines = testPDB.split("\n")
	
	outlines = []
	
	for line in lines:
		if line[:6] in ["ATOM  ", "HETATM", "TER   "]:
			
			l= parseATOMline(line)
			
			if l["ltype"] == "HETATM" and l["rname"] in AA_3_to_1:
				print "HETATM looks like an amino acid:", line,
				l["ltype"] = "ATOM  "
				if l["rname"] == "MET" and l["aname"] == "S":
					l["aname"] = "SD"
					l["aname_exact"] = " SD "
			
			if l["chain"]==" ":
				l["chain"] = "A"
			if l["altloc"]!= " ":
				if l["altloc"] == "A":
					l["altloc"] = " "
				else:
					continue
			if l["rname"] == "MSE":
				l["rname"] = "MET"
				if l["aname"] == "Se":
					l["anme"] = "S"
			
			if l["insert"] != " ":
				continue
			outlines.append(writeATOMdict2Line(l))
		elif line[:6] == "ANISOU":
			continue
		else:
			outlines.append(line)
			
	if not TEST:
		outpdb.write("".join(outlines))
		
		outpdb.close()
	else:
		print "".join(outlines)
	return outpdb.name
	
if __name__=="__main__":
	print run_cleanPDB(sys.argv[1:])
	