import numpy
import math
import random
from io_system import *
from atoms import *
from cell import *

class System(object):
   """
   Represents a simulation cell. 
   Includes the cell parameters, the atoms and the like. """

#__qp holds all the positions and momenta for all the atoms in the simulation
#q and p hold the positions and momenta, respectively.
#P_ext will be the external load.
#The initialisation step now takes a pdc-formatted file for the unit cell and atom positions
#step will eventually call the forces from the external program and then do the propagation step. At the moment we simply take free particle trajectories, to test the theory.
    
   @classmethod
   def from_pdbfile(cls, filedesc, temp = 1.0):
      atoms, cell, natoms = read_pdb(filedesc)
      cls.natoms = natoms
      cls.temp = temp
      cls.k_Boltz = 1.0

      cls.__qp=numpy.zeros((3*natoms,2),float) 
      for i in range(natoms):
         cls.__qp[3*i:3*(i+1),0]=atoms[i][1]
      cls.q=cls.__qp[:,0]

      cls.atoms = [ Atom(cls.__qp[3*i:3*(i+1),:], name = atoms[i][0]) for i in range(natoms) ] #Creates a list of atoms from the __qp array

      cls.P_ext = numpy.zeros((3,3),float)
      cls.cell = Cell(cell, cls.P_ext)

      random.seed(12)
      #cls.__qp[:,1]=numpy.arange(0,3*natoms)*0.01
      for i in range(natoms):
         sigma = math.sqrt(cls.atoms[i].mass * cls.k_Boltz * cls.temp)
         cls.__qp[3*i,1] = random.gauss(0.0, sigma)
         cls.__qp[3*i+1,1] = random.gauss(0.0, sigma)
         cls.__qp[3*i+2,1] = random.gauss(0.0, sigma)
      cls.p=cls.__qp[:,1]
      return cls()

   @classmethod
   def from_system(cls, syst):
      cls.natoms = syst.natoms
      cls.temp = syst.temp
      cls.k_Boltz = syst.k_Boltz

      cls.q = syst.q
      cls.p = syst.p
      cls.__qp = numpy.zeros((3*cls.natoms,2),float) 
      cls.__qp[:,0]=cls.q
      cls.__qp[:,1]=cls.p

      cls.atoms = syst.atoms

      cls.P_ext = syst.P_ext
      cls.cell = syst.cell
      return cls()

   def __str__(self):
      rstr="ATOMS ("+str(self.natoms)+"):\n"
      for i in range(0,self.natoms): 
         rstr=rstr+str(self.atoms[i])+"\n"
      rstr = rstr + "Cell:\n" + str(self.cell)
      rstr = rstr + "\nTotal energy = " + str(self.tot_E()) + ", potential energy = " + str(self.pot()) + ", kinetic energy = " + str(self.kinetic())
      return rstr
       
   def pot(self):
      """Calculates the total potential energy of the system, including
         cell strain"""

      pe = 0.0
      for i in range(self.natoms):
         pe += self.atoms[i].pot()
      pe += self.cell.pot()
      return pe

   def kinetic(self):
      """Calculates the total kinetic energy of the system, including cell 
         kinetic energy"""

      ke = 0.0
      for i in range(self.natoms):
         ke += self.atoms[i].kinetic()
      ke += self.cell.kinetic()
      return ke

   def tot_E(self):
      """Calculates the total energy of the system"""

      return self.kinetic() + self.pot()

#   def step(self,dt):
#      """Takes the atom positions, velocities and forces and integrates the 
#         equations of motion forward by a step dt"""
#      self.q+=self.p*dt

#   def apply_pbc(self):
#      """Takes the system and applies periodic boundary conditions to fold the
#         particle positions back into the unit cell"""
#
#      for i in range(self.natoms):
#         self.atoms[i].q = self.cell.apply_pbc(self.atoms[i])
