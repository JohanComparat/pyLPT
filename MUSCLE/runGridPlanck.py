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
filename="Planck-slice.pdf"
time.time()
p=muscle.generate( ng=256,boxsize=256.)
time.time()
muscle.plotslice(p,filename,boxsize=32.)
f = open('/home2/jcomparat/LPTmeshes/Planck-muscke-2lpt-test.pkl','w')
cPickle.dump(p,f)
f.close()