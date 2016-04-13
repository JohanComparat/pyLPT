import matplotlib
matplotlib.use('pdf')
import pylab as M, numpy as N

def write(redshift=127.,boxsize=512.,ngrid=128,
          hubble=73., ombh2 = 0.045*0.73**2, omch2 = (0.25-0.045)*0.73**2, omk=0.,
          kmax=30.,path='/Users/neyrinck/spikesims/',fileroot='test',
          TimeBetSnapshot=N.sqrt(2.),TimeMax=1.,TimeOfFirstSnapshot=1./16.,Softening=None,SofteningPhys=None):

    if Softening == None:
        Softening = boxsize/ngrid/35.
    if SofteningPhys == None:
        SofteningPhys = boxsize/ngrid/35.
    
    F=open(path+fileroot+'.param','w')
    F.write('InitCondFile  	'+fileroot+'.ic\n')
    F.write('OutputDir          '+fileroot+'.OUT/\n')
    F.write('EnergyFile         energy.txt\n')
    F.write('InfoFile           info.txt\n')
    F.write('TimingsFile        timings.txt\n')
    F.write('CpuFile            cpu.txt\n')

    F.write('RestartFile        restart\n')
    F.write('SnapshotFileBase   '+fileroot+'\n')

    F.write('OutputListFilename outputs.txt\n')

    F.write('% CPU time -limit\n')

    F.write('TimeLimitCPU      345600  % = 96 hours\n')
    F.write('ResubmitOn        0\n')
    F.write('ResubmitCommand   my-scriptfile  \n')
    F.write('\n')
    F.write('\n')
    F.write('% Code options\n')
    F.write('\n')
    F.write('\n')
    F.write('ICFormat                 1\n')
    F.write('SnapFormat               1\n')
    F.write('ComovingIntegrationOn    1\n')
    F.write('\n')
    F.write('TypeOfTimestepCriterion  0\n')
    F.write('OutputListOn             0\n')
    F.write('PeriodicBoundariesOn     1\n')
    F.write('\n')
    F.write('%  Characteristics of run\n')
    F.write('\n')
    F.write('TimeBegin           %g\n'%(1./(1.+redshift)))
    F.write('TimeMax             %g\n'%TimeMax)
    F.write('\n')
    littleh = hubble/100.
    Omega0 = (omch2+ombh2)/littleh**2
    print 'littleh, Omega0 = ',littleh,Omega0,omch2,ombh2
    F.write('Omega0                %g\n'%Omega0)
    F.write('OmegaLambda           %g\n'%(1.-omk-Omega0))
    F.write('OmegaBaryon           %g\n'%(ombh2/littleh**2))
    F.write('HubbleParam           %g\n'%littleh)
    F.write('BoxSize               %g\n'%boxsize)
    F.write('\n')
    F.write('% Output frequency\n')
    F.write('TimeBetSnapshot       '+str(TimeBetSnapshot)+'\n')
    F.write('TimeOfFirstSnapshot   '+str(TimeOfFirstSnapshot)+'\n')
    F.write('\n')
    F.write('CpuTimeBetRestartFile     36000.0    ; here in seconds\n')
    F.write('TimeBetStatistics         0.05\n')
    F.write('\n')
    F.write('NumFilesPerSnapshot       1\n')
    F.write('NumFilesWrittenInParallel 1\n')
    F.write('\n')
    F.write('\n')
    F.write('\n')
    F.write('% Accuracy of time integration\n')
    F.write('\n')
    F.write('ErrTolIntAccuracy      0.025 \n')
    F.write('\n')
    F.write('MaxRMSDisplacementFac  0.2\n')
    F.write('\n')
    F.write('CourantFac             0.15     \n')
    F.write('\n')
    F.write('MaxSizeTimestep       0.03\n')
    F.write('MinSizeTimestep       0.0\n')
    F.write('\n')
    F.write('\n')
    F.write('\n')
    F.write('\n')
    F.write('% Tree algorithm, force accuracy, domain update frequency\n')
    F.write('\n')
    F.write('ErrTolTheta            0.5            \n')
    F.write('TypeOfOpeningCriterion 1\n')
    F.write('ErrTolForceAcc         0.005\n')
    F.write('\n')
    F.write('\n')
    F.write('TreeDomainUpdateFrequency    0.1\n')
    F.write('\n')
    F.write('\n')
    F.write('%  Further parameters of SPH\n')
    F.write('\n')
    F.write('DesNumNgb              60\n')
    F.write('MaxNumNgbDeviation     2\n')
    F.write('ArtBulkViscConst       0.8\n')
    F.write('InitGasTemp            0        % always ignored if set to 0 \n')
    F.write('MinGasTemp             50.0    \n')
    F.write('\n')
    F.write('\n')
    F.write('% Memory allocation\n')
    F.write('\n')
    F.write('PartAllocFactor       1.25\n')
    F.write('TreeAllocFactor       0.8\n')
    F.write('BufferSize            100          % in MByte\n')
    F.write('\n')
    F.write('\n')
    F.write('% System of units\n')
    F.write('\n')
    F.write('UnitLength_in_cm         3.085678e24        ;  1.0 Mpc \n')
    F.write('UnitMass_in_g            1.989e43           ;  1.0e10 solar masses \n')
    F.write('UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec \n')
    F.write('GravityConstantInternal  0\n')
    F.write(' \n')
    F.write('\n')
    F.write('% Softening lengths\n')
    F.write('\n')
    F.write('MinGasHsmlFractional 0.25\n')
    F.write('\n')
    F.write('SofteningGas       0\n')
    F.write('SofteningHalo      '+str(Softening)+'\n')
    F.write('SofteningDisk      0\n')
    F.write('SofteningBulge     0           \n')
    F.write('SofteningStars     0\n')
    F.write('SofteningBndry     0\n')
    F.write('\n')
    F.write('SofteningGasMaxPhys       0.0\n')
    F.write('SofteningHaloMaxPhys      '+str(SofteningPhys)+'\n')
    F.write('SofteningDiskMaxPhys      0\n')
    F.write('SofteningBulgeMaxPhys     0           \n')
    F.write('SofteningStarsMaxPhys     0\n')
    F.write('SofteningBndryMaxPhys     0\n')

def writegadget(pos,vel,redshift,boxsize,OmegaM,OmegaL,HubbleParam, id=None, filename='test.gad',debug=False):

    if len(pos.shape) == 4:
        npart = N.prod(pos.shape[0:3])
    elif len(pos.shape) == 2:
        npart = pos.shape[0]
    
    npartarr = N.array([0,npart,0,0,0,0]).astype(N.int32)
    if debug:
        print pos.shape
        print boxsize,npart
        print 'mass array = ',mass
    mass = N.array([0.,27.757*(boxsize**3/float(npart))*OmegaM,0.,0.,0.,0.]).astype(N.float64)
    time = 1./(1.+redshift) #scale factor

    F= open(filename,'wb')
    F.write(N.array(256).astype(N.int32)) #f77dummy
    F.write(npartarr) # 6-member array
    F.write(mass) # 6-member array
    F.write(N.array([time]))
    F.write(N.array([redshift]))
    F.write(N.array([0,0]).astype(N.int32)) #flagsfr, flag_feedback
    F.write(npartarr)
    F.write(N.array([0,1]).astype(N.int32)) #FlagCooling,NumFiles
    F.write(N.array([boxsize,OmegaM,OmegaL,HubbleParam]))
    bytesleft = 256-6*4 - 6*8 - 8 - 8 - 2*4-6*4 -4 -4 -4*8
    la = N.zeros(bytesleft/4,dtype=N.int32)
    F.write(la)
    F.write(N.array(256).astype(N.int32)) #f77dummy
    
    #pos, vel, id
    F.write(N.array([12*npart]).astype(N.int32)) #f77dummy
    if len(pos.shape) == 4:
        F.write(pos.reshape(N.prod(pos.shape[0:3]),3)) # might have to transpose
    else:
        F.write(pos) # might have to transpose
    F.write(N.array([12*npart]).astype(N.int32)) #f77dummy
    F.write(N.array([12*npart]).astype(N.int32)) #f77dummy
    if len(vel.shape) == 4:
        F.write(vel.reshape(N.prod(vel.shape[0:3]),3)/N.sqrt(time)) # might have to transvele
        # note 1/sqrt(a)
    else:
        F.write(vel) # might have to transpose
    F.write(N.array([12*npart]).astype(N.int32)) #f77dummy

    if id == None:
        id = 1 + N.arange(npart,dtype=N.int32) # Note infuriating "1"!

    F.write(N.array([4*npart]).astype(N.int32)) #f77dummy
    F.write(id) #id array
    F.write(N.array([4*npart]).astype(N.int32)) #f77dummy

    F.close()

    return

def vobozout(pos,filename,f77=False):
    F = open(filename,'w')

    
    np = pos.shape[0]
    print 'np=',np
    N.array([np],dtype=N.int32).tofile(F)
    for d in range(pos.shape[1]):
        (pos[:,d]).astype(N.float32).tofile(F)

    F.close()

