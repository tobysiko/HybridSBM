import sys,numpy
import matplotlib.pyplot as p
import matplotlib.colors as c

def matrixFromFile(f):
	m = []
	for line in open(f).readlines():
		tmp = line.strip().split(",")
		#m.append([float(i) if i != "None" else None for i in tmp ])
		m.append(map(float,line.strip().split(",")))
	return numpy.array(m)






fname1 = sys.argv[1]
fname2 = sys.argv[2]
data1 = matrixFromFile(fname1)
data2 = matrixFromFile(fname2)

maxF=6.0

#diff = data1-data2
diff = [[None for j in xrange(len(data1[i]))]  for i in xrange(len(data1))]
for i in xrange(len(data1)):
	for j in xrange(len(data1[i])):
		if data1[i][j] != maxF and data2[i][j] != maxF:
			diff[i][j] = data1[i][j]-data2[i][j]
		else:
			diff[i][j] = 0.0
			

minimum = numpy.min(diff)
maximum = numpy.max(diff)
print minimum, maximum


ext = 4   #max(abs(minimum),abs(maximum))


print minimum, maximum, ext

ncontours = 400#int(ext*10)+1
step = float(2*ext)/(ncontours)
levels = numpy.arange(-ext,ext+step,step)
#print levels,len(levels),step
p.contourf(diff,cmap='bwr',levels=levels)
#p.contourf(diff,ncontours,cmap='bwr')
p.colorbar(extend="both",ticks=range(-ext,ext+1))
#p.colorbar()
ax = p.gca()
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
p.draw()
p.savefig("diffmap_%s_%s.png"%(fname1[:-4],fname2[:-4]), dpi=400)
p.show()