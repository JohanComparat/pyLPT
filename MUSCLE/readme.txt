Thank you for downloading the MUSCLE (MUltiscale Spherical ColLapse Evolution)
package. It uses updated code taken from CosmoPy (http://www.ifa.hawaii.edu/cosmopy/).
See the MUSCLE paper at http://arxiv.org/abs/1503.07534

Requirements:
1. numpy, matplotlib, scipy.

2. CAMB (available at http://camb.info). Needs to be runnable from 
   ../CAMB/ (relative to running directory of the MUSCLE package). 
   In the current (Feb 2015) version, the CAMB file 
   HighLExtrapTemplate_lenspotentialCls.dat needs to be accessible from the
   local directory, which can be accomplished by changing the path, or
   putting in a link to the file in ../CAMB

MUSCLE has many capabilities: it can generate particle realizations using:
- the Zel'dovich approx
- 2LPT,
- Spherical Collapse (SC) applied only on the interparticle scale, and
- MUltiscale Spherical ColLapse Evolution MUSCLE (default). 
- Scale interpolation with a Gaussian of a specified scale (as in ALPT),
  between any two of these.

A very simple example is in

muscle.example(). This generates a set of 64^3 particles in a box of size
 64 Mpc/h, using MUSCLE.

MUSCLE also can generate simulation initial conditions to be fed into Gadget,
 using muscle.makeic() 
 default method: 2LPT; can also use Zel'dovich
 Other methods currently not supported, since they lack analytic prescriptions for velocities. Differencing nearby time slices could work, though]
