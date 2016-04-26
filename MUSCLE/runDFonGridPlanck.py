import cPickle
import numpy as n
import astropy.cosmology as co
import astropy.units as uu
aa =co.Planck13
import muscle	
import pt
import param
import time
from astropy.io import fits

def prodBox(ng=254, boxsize=64):
	name = "Planck-ng"+str(ng)+"-L"+str(boxsize)
	t0 = time.time()
	p32=muscle.generate(ng = ng, boxsize = boxsize )
	t1 = time.time()
	print "time needed = ", t1 - t0  # 206 s
	muscle.plotslice(p32, name+".pdf",boxsize = boxsize)
	f = open('/home2/jcomparat/LPTmeshes/'+name+'.pkl','w')
	cPickle.dump(p32,f)
	f.close()
	return t1 - t0

ti = time.time()
print ti
z = 0.0
ng=272
boxsize=128.
name = "Planck-ng"+str(ng)+"-L"+str(boxsize)
tracer = 'lrg'
path_to_outputCat = '/home2/jcomparat/LPTmeshes/'+name+'-'+tracer+'.fits.gz'

prodBox(ng, boxsize)

# global deduced parameters 
Lbox = boxsize * uu.megaparsec
massInBox = (aa.critical_density(z) * aa.Om(z) * (Lbox.to(uu.cm))**3).to(uu.solMass)
massPerParticle = massInBox / ng**3.
Mpart = massPerParticle.value
DFunit = (aa.critical_density(z) * aa.Om(z)).to(uu.solMass/uu.megaparsec**3)

# mesh parameters :
dx = boxsize / (ng) # 1000/2048.
cellVolume = (dx*uu.megaparsec)**3

# conversion N particle to DF value
conversion = ((massPerParticle / cellVolume) / (DFunit)).value

# opens the filewith particles
f = open('/home2/jcomparat/LPTmeshes/'+name+'.pkl','r')
p32= cPickle.load(f)
f.close()
p33 = p32.reshape(ng**3,3)

print time.time()
# creation of the mesh
meshsize = int(boxsize / dx) 
DF = n.zeros((meshsize, meshsize, meshsize)) 

# loops over particles and adds them to the mesh
for ii in range(len(p33)):
	index= (p33[ii]/dx).astype(int)
	if len((index==meshsize).nonzero()[0])>0:
		index[(index==meshsize)] = n.ones_like(index[(index==meshsize)])*(meshsize-1)
		DF[index[0],index[1],index[2]]+=1
	else :
		DF[index[0],index[1],index[2]]+=1

DF *= conversion

print time.time()

# now opens the probability distribution of galaxies
px = [8.281210709565847055e-04,
-1.186529082970092658e-02,
6.764603609027315667e-02,
-1.871277409880748865e-01,
2.218162243710562198e-01,
4.102182885734292905e-02, 
-3.649264833679123354e-01, 
3.766102780127327243e-01, 
-4.126106574857090759e-01, 
2.911446087618165812e-01, 
2.280093715491228412e+00, 
-7.436849283207645378e+00]

# random assignment within the cell

catalog=[]
Rs = n.random.uniform(0,1,(meshsize, meshsize, meshsize))
I, J, K = n.meshgrid(n.arange(meshsize), n.arange(meshsize), n.arange(meshsize))
sel = (DF>threshold)

proba = 10**n.polyval(px, n.log10(DF[sel]))
Rs2 = n.random.uniform(0,1,len(proba))
#lrg = (Rs[sel]<proba*10)
lrg = (Rs2<proba)

indices = n.transpose([ I[sel][lrg], J[sel][lrg], K[sel][lrg] ])

minX = I[sel][lrg] * dx
maxX = (I[sel][lrg]+1) * dx
minY = J[sel][lrg] * dx
maxY = (J[sel][lrg]+1) * dx
minZ = K[sel][lrg] * dx
maxZ = (K[sel][lrg]+1) * dx

# from i,j,k to x,y,z
catalog = n.empty((len(minX),4))
for ll in range(len(minX)):
	catalog[ll] = n.random.uniform(minX[ll],maxX[ll]), n.random.uniform(minY[ll],maxY[ll]), n.random.uniform(minZ[ll],maxZ[ll]), DF[sel][lrg][ll]
	
	
"""
threshold = 2.
for ii in range(len(p33)):
	index= (p33[ii]/dx).astype(int)
	if len((index==meshsize).nonzero()[0])>0:
		index[(index==meshsize)] = n.ones_like(index[(index==meshsize)])*(meshsize-1)
		DFval = DF[index[0],index[1],index[2]]
		if DFval < threshold:
			pass
	else :
		DFval = DF[index[0],index[1],index[2]]
		if DFval < threshold:
			pass
	if DFval >=index threshold and Rs[ii]<n.polyval(px, DFVal) :
		catalog.append(n.hstack((p33[ii], DFval)))



catalog=[]
Rs = n.random.uniform(0,1,len(p33))
threshold = 2.
for ii in range(len(p33)):
	index= (p33[ii]/dx).astype(int)
	if len((index==meshsize).nonzero()[0])>0:
		index[(index==meshsize)] = n.ones_like(index[(index==meshsize)])*(meshsize-1)
		DFval = DF[index[0],index[1],index[2]]
		if DFval < threshold:
			pass
	else :
		DFval = DF[index[0],index[1],index[2]]
		if DFval < threshold:
			pass
	if DFval >=index threshold and Rs[ii]<n.polyval(px, DFVal) :
		catalog.append(n.hstack((p33[ii], DFval)))

"""
out = n.transpose(catalog)
		
c0 = fits.Column(name="DF",format='D', array=out[-1] )
c1 = fits.Column(name="x",format='D', array=out[0])
c2 = fits.Column(name="y",format='D', array=out[1] )
c3 = fits.Column(name="z",format='D', array=out[2] )
# now writes the catalog
cols = fits.ColDefs([c1, c2, c3, c4 ])
hdu = fits.BinTableHDU.from_columns(cols)
os.system("rm -rf "+path_to_outputCat)
hdu.writeto(path_to_outputCat)

tf = time.time()

print tf - ti, "seconds, ", (tf - ti)/60., "minutes"