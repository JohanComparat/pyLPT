import matplotlib
matplotlib.use('pdf')
import pylab as M
import numpy as N
import pt
import scipy.signal as SS
import gadgetparam
import os
import itertools

def getgrowth(c,z=0.):
    PS = pt.PowerSpectrum(c.cp)
    return PS.d1(z)

def getkgrid(ng=128,boxsize=512.,whattoget = 'k'):
    thirdim=ng/2+1
    sk = (ng,ng,thirdim)

    kx = N.fromfunction(lambda x,y,z:x, sk).astype(N.float32)
    kx[N.where(kx > ng/2)] -= ng
    kz = N.fromfunction(lambda x,y,z:z, sk).astype(N.float32)
    kz[N.where(kz > ng/2)] -= ng

    if whattoget != 'kxkzk2':
        kmin = 2.*N.pi/boxsize
        kx *= kmin
        kz *= kmin
        
    k2 = kx**2+kx.swapaxes(0,1)**2+kz**2
    if whattoget == 'kxkzOk2':
        k2[0,0,0] = 1.
        return kx/k2, kz/k2
    elif whattoget == 'kxkzk2':
        k2[0,0,0]= 1.
        return kx, kz, k2
    elif whattoget == 'k2':
        return k2
    elif whattoget == 'k':
        return N.sqrt(k2)
    else:
        print 'I dunno what to get'

def sc(dzeld,alpha=1.,mean0 = True):
    """ spherical collapse, applied only at the interparticle scale """

    wgt3=N.where(dzeld > -1.5)
    wle3=N.where(dzeld <= -1.5)
    psi = N.empty(dzeld.shape)
    psi[wgt3] = 3.*(N.sqrt(1.+dzeld[wgt3]*2./3.)-1.)
    psi[wle3] = -3.
    if mean0:
        psi[wgt3] -= N.sum(psi.flatten())/len(wgt3[0])
    return psi

def decreaseres(d,f=2,sumormean='mean'):

    if (f == 1):
        return(d)

    if (len(d.shape) == 3):
        ng = len(d[:,0,0])
        newng = ng/f
        newd = N.zeros(newng**3,dtype=float).reshape((newng,newng,newng))
        for i in range(newng):
            for j in range(newng):
                for k in range(newng):
                    newd[i,j,k] = N.mean(d[f*i:f*(i+1),f*j:f*(j+1),f*k:f*(k+1)])

    elif (len(d.shape) == 2):
        ng = len(d[:,0])
        newng = ng/f
        newd = N.zeros(newng**2.).reshape((newng,newng))
        for i in range(newng):
            for j in range(newng):
                newd[i,j] = N.mean(d[f*i:f*(i+1),f*j:f*(j+1)])

    elif (len(d.shape) == 1):
        ng = len(d)
        newng = ng/f
        newd = N.zeros(newng,dtype=N.float32)
        for i in range(newng):
            newd[i] = N.mean(d[f*i:f*(i+1)])
    else:
        print 'Density array not properly formatted as on a grid'

    return(newd)

def smoothgauss(d,sigma=8.,boxsize=500.,normalize=False,oneminus=False):

    sigmak = 2.*N.pi/sigma
    kmin = 2.*N.pi/boxsize

    dk = N.fft.rfftn(d)
    s = d.shape
    sk = dk.shape
    k2 = getkgrid(s[0],boxsize=boxsize,whattoget='k2')
    if (len(s) == 3):
        a = N.fromfunction(lambda x,y,z:x, sk)
        a[N.where(a > s[0]/2)] -= s[0]
        b = N.fromfunction(lambda x,y,z:y, sk)
        b[N.where(b > s[1]/2)] -= s[1]
        c = N.fromfunction(lambda x,y,z:z, sk)
        c[N.where(c > s[2]/2)] -= s[2]

        k2 = kmin**2*(a**2+b**2+c**2)

    elif (len(s) == 2):
        b = N.fromfunction(lambda y,z:y, sk)
        b[N.where(b > s[0]/2)] -= s[0]
        c = N.fromfunction(lambda y,z:z, sk)
        c[N.where(c > s[1]/2)] -= s[1]

        k2 = kmin**2*(b**2+c**2)
    elif (len(s) == 1):
        c = N.fromfunction(lambda z:z, sk)
        c[N.where(c > s[0]/2)] -= s[0]
        
        k2 = kmin**2*c**2
        
    gaussian = N.exp(-k2/sigmak**2/2.)
    if oneminus:
        gaussian = 1.-gaussian
    ans = N.fft.irfftn(dk*gaussian)
    if normalize:
        ans /= N.sqrt(2.*N.pi*sigmak**2)

    return ans

def muscle(dzeld, alpha=1.,mean0=False,smooth='Gauss'):
    """ MUltiscale Spherical colLapse Evolution """
    ng=dzeld.shape[0]
    twofolds = int(N.log(ng)/N.log(2.))
    comparator = 1+(2./3.)*dzeld

    print 'min(comparator)=',N.min(comparator)
    if smooth == 'degrade':
        starter = 1
    elif smooth == 'Gauss':
        starter = 0
    for i in N.arange(starter,twofolds):
        f=2**i
        if smooth == 'degrade':
            dzeld_degrade = N.repeat(N.repeat(N.repeat(
                        distrib.decreaseres(dzeld,f=f),
                        f,axis=0),f,axis=1),f,axis=2)
            # degrades the resolution by a factor of 2 in each direction,
            # and then produces a grid with the original dimensions, with
            # 2x2x2 supercells

        elif smooth == 'Gauss':
            dzeld_degrade = smoothgauss(dzeld,sigma=f,boxsize=ng)
        comparator_degrade = 1+(2./3.)*dzeld_degrade
        mincomp = N.min(comparator_degrade)
        print 'mincomp=',mincomp
        if mincomp > 0.: # if we're so low-res that nothing's collapsing
            break
        w=N.where((comparator_degrade <= 0.) | (comparator <= 0.))
        comparator[w] = N.minimum(comparator_degrade[w],comparator[w])

    psi = -3.*N.ones(dzeld.shape,dtype=float)
    wherenocollapse = N.where(comparator > 0.)
    psi[wherenocollapse] = 3.*(N.sqrt(1+(2./3.)*dzeld[wherenocollapse])-1.)

    return psi

def makegauss(camby,kgrid,seed=314159,boxsize=None,exactpk=False,returnfft=False):
    """ 
    Makes a gaussian random field from a Camb instance (camby), which
      holds the linear power spectrum
    exactpk -- if True, sets Fourier amplitudes exactly to the linear power spectrum
    returnfft -- if True, returns the FFT of the density field
    """
    ng = kgrid.shape[0]
    thirdim = ng/2+1

    rs = N.random.RandomState(seed)

    if (exactpk):
        dk = N.exp(2j*N.pi*rs.rand(ng*ng*(thirdim))).reshape((ng,ng,thirdim)).astype(N.complex64)
    else:
        dk = N.empty((ng,ng,thirdim),dtype=N.complex64)
        dk.real = rs.normal(size=ng*ng*(thirdim)).reshape((ng,ng,thirdim)).astype(N.float32)
        dk.imag = rs.normal(size=ng*ng*(thirdim)).reshape((ng,ng,thirdim)).astype(N.float32)
        dk /= N.sqrt(2.)

    sk = kgrid.shape
    dk *= N.sqrt(camby.pkInterp(kgrid.flatten())).reshape(sk)*ng**3/boxsize**1.5

    # Nyquist crap: dk(-k) = conjugate(dk(k))
    dk[ng/2+1:,1:,0]= N.conj(N.fliplr(N.flipud(dk[1:ng/2,1:,0])))
    dk[ng/2+1:,0,0] = N.conj(dk[ng/2-1:0:-1,0,0])
    dk[0,ng/2+1:,0] = N.conj(dk[0,ng/2-1:0:-1,0])
    dk[ng/2,ng/2+1:,0] = N.conj(dk[ng/2,ng/2-1:0:-1,0])

    dk[ng/2+1:,1:,ng/2]= N.conj(N.fliplr(N.flipud(dk[1:ng/2,1:,ng/2])))
    dk[ng/2+1:,0,ng/2] = N.conj(dk[ng/2-1:0:-1,0,ng/2])
    dk[0,ng/2+1:,ng/2] = N.conj(dk[0,ng/2-1:0:-1,ng/2])
    dk[ng/2,ng/2+1:,ng/2] = N.conj(dk[ng/2,ng/2-1:0:-1,ng/2])

    if (returnfft):
        return dk
    else:
        return N.fft.irfftn(dk)

def psi2pos(psi,x=None,boxsize=None,dopbc=True,center=True,xdtype=N.float32):
    """ makes a position field from a displacement field """
    ng = psi.shape[0]
    pos = 1.*psi
    x=boxsize/ng*N.fromfunction(lambda x,y,z:x, (ng,ng,ng))
    pos[:,:,:,0] += x
    pos[:,:,:,1] += x.swapaxes(0,1)
    pos[:,:,:,2] += x.swapaxes(0,2)
    if dopbc:
        pos = pos%boxsize
    return pos

def scalegaussinterp(dlowk,dhighk,kgrid,sigma=5.):
    """ 
    Interpolates between large- and small-scale displacement divergences.
    """

    ng = kgrid.shape[0]
    boxsize=2.*N.pi/kgrid[1,0,0]
    print boxsize, sigma
    print 'ratio of boxsize to sigma:',boxsize/sigma
    print 'ratio of sigma to cellsize:',sigma/(boxsize/ng)
    sigmak = 2.*N.pi/sigma
    gaussian = N.exp(-(kgrid/sigmak)**2/2.)

    return dlowk*gaussian+dhighk*(1.-gaussian)

def invdiv(d=None,dk=None,boxsize=None,dtype=N.float32,dopsi2pos=True,dopbc=True):
    """ Fourier-space inverse-divergence """
    # -i\bk/k^2 
    if ((d == None) & (dk == None)):
        print 'we have to have some input!'

    if dk == None:
        dk = N.fft.rfftn(d)
    ng = dk.shape[0]

    kx,kz = getkgrid(ng=ng,boxsize=boxsize,whattoget='kxkzOk2')
    pos = N.zeros((ng,ng,ng,3),dtype=dtype)
    pos[:,:,:,0] = N.fft.irfftn(-1j*dk*kx)
    pos[:,:,:,1] = N.fft.irfftn(-1j*dk*kx.swapaxes(0,1))
    pos[:,:,:,2] = N.fft.irfftn(-1j*dk*kz)

    if dopsi2pos:
        pos = psi2pos(pos,boxsize=boxsize,dopbc=dopbc)
    return pos

def twolpt(dpk,boxsize=None):
    """ Takes Zel'dovich (div_L . psi)_k """
    ng = dpk.shape[0]
    kx,kz,k2 = getkgrid(ng=ng,boxsize=boxsize,whattoget='kxkzk2')

    growth2=-3./7.#*lagrexp.getgrowth(z=z)

    # first the diagonal
    phi0 = N.fft.irfftn(dpk*kx**2/k2)
    phi1 = N.fft.irfftn(dpk*kx.swapaxes(0,1)**2/k2)
    phi2 = N.fft.irfftn(dpk*kz**2/k2)

    ans = phi0*phi1+phi0*phi2+phi1*phi2
    phi0 = None; phi1 = None; phi2 = None

    #phi01, phi02,phi12
    ans -= N.fft.irfftn(dpk*kx*kx.swapaxes(0,1)/k2)**2
    ans -= N.fft.irfftn(dpk*kx*kz/k2)**2
    ans -= N.fft.irfftn(dpk*kx.swapaxes(0,1)*kz/k2)**2

    # Needs dp=ifft(dpk) to be added to it! We do that afterward, having saved dp
    return growth2*ans

def f_omega(c,a):
    Omega_M = (c.cp.ombh2+c.cp.omch2)/(c.cp.hubble/100.)**2
    Omega_L = 1.-c.cp.omk - Omega_M
    Omega_a = Omega_M/(Omega_M + a*c.cp.omk + a**3 * Omega_L)
    return Omega_a**0.6

def f2_omega(c,a):
    Omega_M = (c.cp.ombh2+c.cp.omch2)/(c.cp.hubble/100.)**2
    Omega_L = 1.-c.cp.omk - Omega_M
    Omega_a = Omega_M/(Omega_M + a*c.cp.omk + a**3 * Omega_L)
    return 2.*Omega_a**(4./7.)

def generate(c=None,d=None,dk=None,ng=64,boxsize=128.,sigma8=0.829, sigmaalpt = 10.,scheme='muscle',largescale=None,smallscale=None,exactpk=False, seed = 42389,returnfft = False, dopbc=True, fileroot=None, returnpos=True, returndisp = False, plottest=False, returninitdens=False, justreturnc = False, returnvel=False,deltaz4vel=1./128., hubble=67.77, ombh2 = 0.048252*0.6777**2, omch2 = (0.30712-0.048252)*0.6777**2, redshift = 0.,kmax=30.,omk=0.):
    """ 
    possible inputs:
    c = Camb instance; contains cosmological parameters, power spectrum
    d = configuration-space density field
    dk = FFT of a density field

    parameters:
    sigmaalpt = Gaussian k-space interpolation smoothing scale, as in ALPT
    scheme: can be 'zeld'ovich, '2lpt', 'sc' (single-scale spherical collapse), 'muscle' (multiscale)
    largescale/smallscale: use for scale interpolation
    dopbc: False to preserve distortion of particle lattice, not enforcing periodic boundary conditions
    returnpos, returndisp, returnvel: return position, displacement field at particles, velocities [velocities only work for Zeld, 2LPT currently!]
    plottest: show a slice of the particles
    exactpk: Set each Fourier amplitude exactly to the linear power spectrum, 
             suppressing fluctuations in Fourier amplitudes
    """
    
    if (returnpos & returndisp):
        print 'Choose either position or displacement field to return'
        return
    if returnvel:
        omegam = (ombh2+omch2)/(hubble/100.)**2
        omegal = 1.-omk-omegam
        hubble_z = 100.*N.sqrt(omegam*(1.+redshift)**3 + omk*(1+redshift)**2 + omegal)
        # Valid only for simple dark energy; doesn't allow w != -1

    if ((((c == None) | (returnvel & (smallscale != None))) &
        ((d == None) & (dk == None))) | (justreturnc == True)):
        print 'c == None:', (c == None)
        print 'returnvel = ', returnvel
        c = pt.Camb(hubble=hubble, ombh2 = ombh2, omch2 = omch2,omk=omk,
                    transfer_kmax=kmax,transfer_redshift=[0.])
        if (dk == None):
            c.run()
            sigma82default_z0 = pt.normalizePk(c,sigma8)
            if redshift != 0.:
                c = pt.Camb(hubble=hubble, ombh2 = ombh2, omch2 = omch2,omk=omk,
                            transfer_kmax=kmax,transfer_redshift=[redshift])
                c.run()
                c.pk *= sigma8**2/sigma82default_z0
                c.logpk = N.log(c.pk)
                c.pkSplineCoeff = SS.cspline1d(c.logpk)
        if justreturnc:
            return c

        if (returnvel & (smallscale != None)):
            print 'Numerically computing velocity'
            cplusdeltaz = pt.Camb(hubble=hubble, ombh2 = ombh2, omch2 = omch2, omk=omk,
                                  transfer_kmax=kmax,transfer_redshift=[redshift+deltaz4vel])
            cplusdeltaz.run()

            cplusdeltaz.pk *= sigma8**2/sigma82default_z0
            cplusdeltaz.logpk = N.log(cplusdeltaz.pk)
            growthindeltaz=N.mean(c.logpk-cplusdeltaz.logpk)

            print 'camb:', growthindeltaz
            print 'naive:',2.*N.log(getgrowth(c,z=redshift)/getgrowth(c,z=redshift+deltaz4vel))
        
            noiselevel = N.std(c.logpk-cplusdeltaz.logpk)
            print 'std/mean(growth from camb)=',noiselevel/growthindeltaz
            if (noiselevel/growthindeltaz > 1./8.):
                print "Warning! deltaz so small that it's giving lots of noise."
            cplusdeltaz.pk = c.pk * N.exp(growthindeltaz)
            cplusdeltaz.logpk = c.logpk + growthindeltaz
            cplusdeltaz.pkSplineCoeff = SS.cspline1d(cplusdeltaz.logpk)
            
            #The Hubble parameter at z=redshift
            dispplusdeltaz = generate(c=cplusdeltaz,ng=ng,boxsize=boxsize,sigmaalpt = sigmaalpt,
                                     largescale=largescale,smallscale=smallscale, seed = seed,
                                     dopbc=dopbc, fileroot=fileroot, returnpos=False, returndisp=True,
                                     plottest=False, returnvel=False,exactpk=exactpk)

    # the Zel'dovich displacement-divergence in Fourier space
    if d == None:
        if dk == None:
            kgrid = getkgrid(ng=ng,boxsize=boxsize,whattoget='k')
            dk=makegauss(c,kgrid,boxsize=boxsize,exactpk=exactpk,returnfft=True,seed=seed)
            if returnfft:
                return dk
        else:
            #(over)write ng
            ng=dk.shape[0]
        d = N.fft.irfftn(dk)

    else: #shouldn't have both d and dk non-None
        print 'd supplied'
        ng = d.shape[0]
        kgrid = getkgrid(ng=ng,boxsize=boxsize,whattoget='k')
        dk = N.fft.rfftn(d)

    if ((scheme != None) & ((smallscale != None) | (largescale != None))):
        print "Please specify only 'scheme' or ('smallscale' and 'largescale')"
    
    if (smallscale == None) & (largescale == None):
        if (scheme != None):
            largescale = scheme
        else:
            print "Please specify either a full 'scheme' or "
            print "a 'smallscale' and 'largescale' displacement field scheme."
            return

    print 'largescale=',largescale

    if smallscale == 'sc':
        psismall = N.fft.rfftn(sc(-d))
    elif smallscale == 'muscle':
        psismall = N.fft.rfftn(muscle(-d))
    elif smallscale == 'zeld':
        psismall = -dk

    if largescale == 'zeld':
        psilarge = -dk
    elif largescale == '2lpt':
        psiquadratic = N.fft.rfftn(twolpt(-dk,boxsize=boxsize))
        psilarge = psiquadratic - dk
        # dk because we didn't add it in twolpt for efficiency reasons
    elif largescale == 'sc':
        psilarge = N.fft.rfftn(sc(-d))
    elif largescale == 'muscle':
        psilarge = N.fft.rfftn(muscle(-d))

    if (smallscale != None) & (largescale != None):
        psik = scalegaussinterp(psilarge,psismall,kgrid,sigma=sigmaalpt)
    elif smallscale != None:
        psik = psismall
    elif largescale != None:
        psik = psilarge

    disp = invdiv(dk=psik,boxsize=boxsize,dopsi2pos=False).reshape(ng,ng,ng,3)
    pos = psi2pos(disp,boxsize=boxsize,dopbc=dopbc)

    if returnvel: # only works for Zeld, 2LPT
        time = 1./(1.+redshift)
        print 'time, hubble_z, f_omega = ',time,hubble_z,f_omega(c,time)
        vel = disp * time *hubble_z * f_omega(c,time)
        print 'total factor = ',time*hubble_z*f_omega(c,time)
        print 'factor in gadget = ',N.sqrt(time)*hubble_z*f_omega(c,time)
        if scheme == '2lpt':
            vel += 3./7. * time * hubble_z * f2_omega(c,time)* \
                invdiv(dk=psiquadratic,boxsize=boxsize,dopsi2pos=False)#.swapaxes(0,2)

    if plottest:
        if fileroot == None:
            fileroot='plottest'
        plotslice(pos,filename=fileroot+'.png',boxsize=boxsize)

    if returndisp:
        return disp
    if (returninitdens & (returnpos == False) & (returnvel == False)):
        return d
    if (returnpos & returnvel & (returninitdens == False)):
        return pos,vel
    if (returnpos & returninitdens & (returnvel == False)):
        return pos,d
    if (returnpos & returnvel & returninitdens):
        return pos,vel,d
    if returnpos:
        return pos
    if returnvel:
        return vel
    return

def makeic(c=None,d=None,dk=None, scheme='2lpt', largescale=None, 
           smallscale=None,boxsize=None,ng=None,
           hubble=73., ombh2 = 0.045*0.73**2, omch2 = (0.25-0.045)*0.73**2,sigma8=0.8, 
           kmax=30.,omk=0.,path='./',fileroot='test',deltaz4vel=1./16.,TimeBetSnapshot=2.,pretendEdS=False,
           exactpk=False,returninitdens=False,seed=42389,returnpos=False,redshift=127.,TimeMax=1.,Softening=None,SofteningPhys=None):
    """
    Generate Gadget files suitable for running
    See generate() for argument explanations
    """

    output=generate(c=c,d=d,redshift=redshift,dk=dk,boxsize=boxsize,ng=ng,
                    hubble=hubble, ombh2 = ombh2, omch2 = omch2, omk=omk,
                    kmax=kmax, returnpos=True,returnvel=True,deltaz4vel=deltaz4vel,
                    scheme=scheme,largescale=largescale,smallscale=smallscale,
                    returninitdens=returninitdens,exactpk=exactpk,seed=seed,sigma8=sigma8)

    if returninitdens:
        pos,vel,initdens = output
    else:
        pos,vel = output

    print 'maxpos=',max(pos.flatten())

    os.system('mkdir '+path+fileroot+'.OUT/')
    OmegaM = (ombh2+omch2)/(hubble/100.)**2
    OmegaL = 1.-omk-OmegaM
    gadgetparam.write(boxsize=boxsize,redshift=redshift, hubble=hubble, 
                      ombh2 = ombh2, omch2 = omch2, omk=omk,kmax=kmax,
                      path=path, fileroot=fileroot,TimeBetSnapshot=TimeBetSnapshot, TimeMax=TimeMax)        
    gadgetparam.writegadget(pos,vel,redshift,boxsize,OmegaM,OmegaL,hubble, id=None, filename=path+fileroot+'.ic')

    if returninitdens & returnpos:
        return initdens, pos
    if returninitdens:
        return initdens
    if returnpos:
        return pos

def plotslice(pos,filename='',boxsize=100.):
    ng = pos.shape[0]
    M.clf()
    M.scatter(pos[ng/4,:,:,1].flatten(),pos[ng/4,:,:,2].flatten(),s=1.,lw=0.)
    M.axis('tight')
    if filename != '':
        M.savefig(filename)

def example(filename="example-slice.png"):
    p=generate(ng=64,boxsize=64.)
    plotslice(p,filename)
