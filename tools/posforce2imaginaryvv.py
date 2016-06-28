#!/usr/bin/python
""" posforce2imaginary.py

Reads beads positions and forces from an i-PI run and computes the
imaginary time velocity autocorrelation functions for each particle,
in a single file

Assumes the input files are in xyz format and atomic units,
with prefix.pos_*.xyz and prefix.for_*.xyz naming scheme.

Syntax:
   posforce2imaginaryvv.py prefix temperature[K] 
"""

import numpy as np
import sys, glob
from ipi.utils.io.io_xyz import *
from ipi.engine.beads import Beads
from ipi.utils.depend import *
from ipi.utils.units import *

def main(prefix, temp):

   temp = unit_to_internal("energy","kelvin",float(temp))

   ipos=[]
   for filename in sorted(glob.glob(prefix+".pos*")):
      ipos.append(open(filename,"r"))

   ifor=[]
   for filename in sorted(glob.glob(prefix+".for*")):
      ifor.append(open(filename,"r"))

   ivacfim=open(prefix+".imvacf","w")

   nbeads = len(ipos)

   if (nbeads!=len(ifor)): raise ValueError("Mismatch between number of output files for forces and positions")
   natoms = 0
   ifr = 0
   while True:
      try:
         for i in range(nbeads):
            pos = read_xyz(ipos[i])
            force = read_xyz(ifor[i])
            if natoms == 0:
               natoms = pos.natoms
               beads = Beads(natoms,nbeads)
               forces = Beads(natoms,nbeads)
            beads[i].q = pos.q
            beads[i].m = pos.m
            forces[i].q = force.q
      except EOFError: # finished reading files
         break

      if (ifr==0):
         cim = np.zeros((natoms,nbeads+1),float)

      q = depstrip(beads.q)
      f = depstrip(forces.q)
      m = depstrip(beads.m)
      qc = depstrip(beads.qc)

      for j in range(nbeads+1):
         for i in range(natoms):
#            blablis=[((j+k)%nbeads, k) for k in range(nbeads)]
#            print j, blablis
            dummy=[(q[(j+k)%nbeads, 3*i:3*(i+1)]-qc[3*i:3*(i+1)])*(-f[k, 3*i:3*(i+1)]) for k in range(nbeads)]
            cim[i,j]=cim[i,j]+np.sum(dummy)/(nbeads*m[i])
            cim[i,j]+=3.*Constants.kb*temp/m[i] 

      ifr+=1
      print 'Frame: ', ifr

   cim[:,:]=cim[:,:]/ifr
   
   for j in range(nbeads+1):
      fullcim[j]=np.sum([i for i in cim[:,j]]) 

   ivacfim.write("# Imaginary time autocorrelation, averaged through the whole simulation -- beta*hbar is %f \n" % (Constants.hbar/(Constants.kb*temp)))

   for islice in range(nbeads+1):
      ivacfim.write("%f %12.5e \n" % (float(islice)/float(nbeads), fullcim[islice]))

#   for iatom in range(natoms):
#      ivacfim.write("\n# Imvacf for atom %d %s \n" % (iatom, pos.names[iatom]))
#      for islice in range(nbeads+1):
#         ivacfim.write("%f %12.5e \n" % (float(islice)/float(nbeads), cim[iatom, islice]))

# could add the full imvacf in another file, as well as other types of imvacf 
      

if __name__ == '__main__':
   main(*sys.argv[1:])
