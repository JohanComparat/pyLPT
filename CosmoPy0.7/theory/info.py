""">>> info.py <<<

This module contains functions used to calculate Fisher information-content
curves from power spectrum covariance

Current revision:
    ID:         $Id: info.py 42 2005-11-11 02:49:17Z neyrinck $
    Date:       $Date: 2005-11-10 16:49:17 -1000 (Thu, 10 Nov 2005) $
    Revision:   $Revision: 42 $

(C) 2006 The CosmoPy Team (see Copyright for details)
"""

import pylab as M

def getParamCovMat(prefix,dlogpower = 2, theoconstmult = 1.,dlogfilenames = ['dlogpnldloga.dat'],volume=256.**3,startki = 0, endki = 0, veff = [0.]):
    """
    Calculates parameter covariance matrix from the power spectrum covariance matrix and derivative term
    in the prefix directory
    """
    nparams = len(dlogfilenames)

    kpnl = M.load(prefix+'pnl.dat')
    k = kpnl[startki:,0]

    nk = len(k)
    if (endki == 0):
        endki = nk
        
    pnl = M.array(kpnl[startki:,1],M.Float64)
    covarwhole = M.load(prefix+'covar.dat')
    covar = covarwhole[startki:,startki:]
    if len(veff) > 1:
        sqrt_veff = M.sqrt(veff[startki:])
    else:
        sqrt_veff = M.sqrt(volume*M.ones(nk))

    dlogs = M.reshape(M.ones(nparams*nk,M.Float64),(nparams,nk))
    paramFishMat = M.reshape(M.zeros(nparams*nparams*(endki-startki),M.Float64),(nparams,nparams,endki-startki))
    paramCovMat = paramFishMat * 0.

    # Covariance matrices of dlog's
    for param in range(nparams):
        if len(dlogfilenames[param]) > 0:
            dlogs[param,:] = M.load(prefix+dlogfilenames[param])[startki:,1]

    normcovar = M.zeros(M.shape(covar),M.Float64)
    for i in range(nk):
        normcovar[i,:] = covar[i,:]/(pnl*pnl[i])

    M.save(prefix+'normcovar.dat',normcovar)

    f = k[1]/k[0]

    if (volume == -1.):
        volume = (M.pi/k[0])**3

    #theoconst = volume * k[1]**3 * f**(-1.5)/(12.*M.pi**2) #1 not 0 since we're starting at 1
    for ki in range(1,endki-startki):
        for p1 in range(nparams):
            for p2 in range(nparams):
                paramFishMat[p1,p2,ki] = M.sum(M.sum(\
                M.inverse(normcovar[:ki+1,:ki+1]) *
                M.outerproduct(dlogs[p1,:ki+1]*sqrt_veff[:ki+1],\
                               dlogs[p2,:ki+1]*sqrt_veff[:ki+1])))
                
                
        paramCovMat[:,:,ki] = M.inverse(paramFishMat[:,:,ki])

    return k[1:],paramCovMat[:,:,1:]
