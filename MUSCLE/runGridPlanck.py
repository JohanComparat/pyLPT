"""
c=None,d=None,dk=None,ng=64,boxsize=128.,sigma8=0.829, sigmaalpt = 10.,scheme='muscle',largescale=None,smallscale=None,exactpk=False, seed = 42389,returnfft = False, dopbc=True, fileroot=None, returnpos=True, returndisp = False, plottest=False, returninitdens=False, justreturnc = False, returnvel=False,deltaz4vel=1./128., hubble=67.77, ombh2 = 0.048252*0.6777**2, omch2 = (0.30712-0.048252)*0.6777**2, redshift = 0.,kmax=30.,omk=0.):
possible inputs:
c = Camb instance; contains cosmological parameters, power spectrum
d = configuration-space density field
dk = FFT of a density field

parameters:
sigmaalpt = Gaussian k-space interpolation smoothing scale, as in ALPT
scheme: can be 'zeld'ovich, '2lpt', 'sc' (single-scale spherical collapse), 'muscle' (multiscale)
largescale/smallscale: use for scale interpolationv
dopbc: False to preserve distortion of particle lattice, not enforcing periodic boundary conditions
returnpos, returndisp, returnvel: return position, displacement field at particles, velocities [velocities only work for Zeld, 2LPT currently!]
plottest: show a slice of the particles
exactpk: Set each Fourier amplitude exactly to the linear power spectrum, 
suppressing fluctuations in Fourier amplitudes
"""
import cPickle
import muscle	
import pt
import param
import time
# c = pt.Camb(cambParam = param)

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

prodBox(ng=490, boxsize=128.) # 809 : MDPL R
prodBox(ng=254, boxsize=64.) # 185 : MDPL R

prodBox(ng=328, boxsize=128.) # - : MDPL R
prodBox(ng=164, boxsize=64.) # 26 : 2/3 MDPL R

prodBox(ng=246, boxsize=128.) # 88 : half MDPL R
prodBox(ng=122, boxsize=64.) # 12 : half MDPL R
	

t0 = time.time()
p=muscle.generate( ng=256,boxsize=256.)
t1 = time.time()
print "time needed = ", t1 - t0 # 124 s

t0 = time.time()
p=muscle.generate( ng=512,boxsize=512.)
t1 = time.time()
print "time needed = ", t1 - t0 # 1603.97744703

# p contains the final position of each initial grid point

# then on eneeds to convert to an overdensity and add the biasing model

muscle.plotslice(p,filename,boxsize=32.)
f = open('/home2/jcomparat/LPTmeshes/Planck-muscke-2lpt-test.pkl','w')
cPickle.dump(p,f)
f.close()
