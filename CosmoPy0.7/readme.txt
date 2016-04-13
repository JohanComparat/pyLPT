Hello, and welcome to CosmoPy v. 0.6.  It was written by Mark
Neyrinck, Istvan Szapudi, Adrian Pope, Peter Papai, and Tamas
Budavari.

CosmoPy requires Python, a few packages for it, and CAMB.

The Python packages you will need, in decreasing order of importance
(no pun intended), are:

numpy http://numpy.scipy.org/
pylab (matplotlib) http://matplotlib.sourceforge.net/
scipy http://www.scipy.org/
more optionally:
pygsl http://pygsl.sourceforge.net/

along with these packages which are probably included with Python,
again in rough decreasing order of importance:

copy, unittest, cex, RandomArray

To start using CosmoPy (the following assumed you unzipped things into
the directory ~):

- Make sure Python and the packages mentioned above are installed.

- Download and compile CAMB (http://www.camb.info).  Put a copy of the
CAMB executable ("camb") in ~/CosmoPy0.x/CAMB .  CosmoPy has been
thoroughly tested on the June 2006 version of CAMB; the Sep 2006
version should be fully compatible, too.

- go into ~/CosmoPy0.x/theory, and start python.  Here is a sample
  transcript:

 ~/CosmoPy0.6/theory $ python
 Python 2.4.4 (#1, Jun 15 2007, 13:45:10) 
 [GCC 4.1.1 (Gentoo 4.1.1-r1)] on linux2
 Type "help", "copyright", "credits" or "license" for more information.
 >>> import example
 >>> example.hod()
 Age of universe/GYr  =  13.461
 Om_b h^2             = 0.02450
 Om_c h^2             = 0.12250
 Om_nu h^2            = 0.00000
 Om_Lambda            = 0.70000
 Om_K                 = 0.00000
 Sharp opt depth      =  0.101
 tau_recomb/Mpc       =  278.83  tau_now/Mpc =  13894.7
  at z =   0.0000000E+00  sigma8 (all matter)=  0.8857116    
 delta =  1.67454254708
 Normalization const (integrated): 0.322130600942
 >>>

 If everything is set up correctly, you should get two sets of plots
 of galaxy power spectra for different HOD's.  Also try the other
 example programs, example.redshift() and example.getInfoCurve().

If you need assistance, or if you find a bug, please go to
http://www.ifa.hawaii.edu/cosmopy/.  Soon, we hope to have a
bug-tracking system and wiki set up there.  If that's not set up yet,
(or even if it is) feel free to e-mail Istvan Szapudi
<szapudi@ifa.hawaii.edu> or Mark Neyrinck <neyrinck@ifa.hawaii.edu>.

This document last updated Feb 22, 2008.
