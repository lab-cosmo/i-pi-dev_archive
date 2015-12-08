"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io.backends.io_xyz import read_xyz
from ipi.utils.io.backends.io_pdb import read_pdb
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import Constants, unit_to_internal
from ipi.inputs.thermostats import InputThermo
from ipi.inputs.barostats import InputBaro
from ipi.engine.mover import Mover
from ipi.engine.thermostats import *
from ipi.engine.barostats import *


__all__ = ['DynMover', 'NVEEnsemble', 'NVTEnsemble', 'NPTEnsemble', 'NSTEnsemble']


class DynMover(Mover):
   """Base (path integral) molecular dynamics class.

   Gives the standard methods and attributes needed in all the
   dynamics classes.

   Attributes:
      beads: A beads object giving the atoms positions.
      cell: A cell object giving the system box.
      forces: A forces object giving the virial and the forces acting on
         each bead.
      prng: A random number generator object.
      nm: An object which does the normal modes transformation.      

   Depend objects:
      econs: The conserved energy quantity appropriate to the given
         ensemble. Depends on the various energy terms which make it up,
         which are different depending on the ensemble.he
      temp: The system temperature.
      dt: The timestep for the algorithms.
      ntemp: The simulation temperature. Will be nbeads times higher than
         the system temperature as PIMD calculations are done at this
         effective classical temperature.
   """

   def __init__(self, timestep, mode="nve", edyn=0.0, thermostat=None, barostat = None, fixcom=False, fixatoms=None):
      """Initialises a "dynamics" mover.

      Args:
         dt: The timestep of the simulation algorithms.
         temp: The temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(DynMover, self).__init__(fixcom = fixcom, fixatoms = fixatoms)
      dset(self, "dt",    depend_value(name='dt',    value=timestep))
      dset(self, "edyn",    depend_value(name='edyn',    value=edyn))
      dget(self, "econs")._func = self.get_econs
      if thermostat == None:
         self.thermostat = Thermostat()
      else:
         self.thermostat = thermostat
      
      if barostat == None:
         self.barostat = Barostat()
      else:
         self.barostat = barostat
         
      self.enstype = mode
      self.fixcom = fixcom
      if fixatoms is None:
         self.fixatoms = np.zeros(0,int)
      else:
         self.fixatoms = fixatoms


   def bind(self, ens, beads, nm, cell, bforce, bbias, prng):
      """Binds ensemble beads, cell, bforce, bbias and prng to the dynamics.

      This takes a beads object, a cell object, a forcefield object and a
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Note that the conserved
      quantity is defined in the init, but as each ensemble has a different
      conserved quantity the dependencies are defined in bind.

      Args:
         beads: The beads object from whcih the bead positions are taken.
         nm: A normal modes object used to do the normal modes transformation.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """

      # store local references to the different bits of the simulation      
      self.beads = beads
      self.cell = cell
      self.forces = bforce
      self.bias = bbias
      self.prng = prng
      self.nm = nm
      self.ensemble = ens

      # n times the temperature
      dset(self,"ntemp", depend_value(name='ntemp',func=self.get_ntemp,
         dependencies=[dget(self.ensemble,"temp")]))

      # dependencies of the conserved quantity

      dget(self,"econs").add_dependency(dget(self.beads, "kin"))
      dget(self,"econs").add_dependency(dget(self.forces, "pot"))
      dget(self,"econs").add_dependency(dget(self.bias, "pot"))
      dget(self,"econs").add_dependency(dget(self.beads, "vpath"))
      dget(self,"econs").add_dependency(dget(self.ensemble, "eens"))
      dget(self,"econs").add_dependency(dget(self, "edyn"))
      self.pconstraints() # applies momentum constraints to initial configurations


   def get_ntemp(self):
      """Returns the PI simulation temperature (P times the physical T)."""

      return self.temp*self.beads.nbeads


   def pstep(self):
      """Dummy momenta propagator which does nothing."""

      pass

   def qcstep(self):
      """Dummy centroid position propagator which does nothing."""

      pass

   def step(self, step=None):
      """Dummy simulation time step which does nothing."""

      pass

   def get_econs(self):
      """Calculates the conserved energy quantity for constant energy
      ensembles.
      """
      eham = self.beads.vpath*self.nm.omegan2 + self.nm.kin + self.forces.pot
      eham += self.bias.pot # bias
      return eham + self.ensemble.eens + self.edyn

   def pconstraints(self):
      pass


#~ 
#~ class NVEEnsemble(Ensemble):
   #~ """Ensemble object for constant energy simulations.
#~ 
   #~ Has the relevant conserved quantity and normal mode propagator for the
   #~ constant energy ensemble. Note that a temperature of some kind must be
   #~ defined so that the spring potential can be calculated.
#~ 
   #~ Attributes:
      #~ ptime: The time taken in updating the velocities.
      #~ qtime: The time taken in updating the positions.
      #~ ttime: The time taken in applying the thermostat steps.
#~ 
   #~ Depend objects:
      #~ econs: Conserved energy quantity. Depends on the bead kinetic and
         #~ potential energy, and the spring potential energy.
   #~ """
#~ 
   #~ def __init__(self, dt, temp, fixcom=False, eens=0.0, fixatoms=None):
      #~ """Initialises NVEEnsemble.
#~ 
      #~ Args:
         #~ dt: The simulation timestep.
         #~ temp: The system temperature.
         #~ fixcom: An optional boolean which decides whether the centre of mass
            #~ motion will be constrained or not. Defaults to False.
      #~ """
#~ 
      #~ super(NVEEnsemble,self).__init__(dt=dt,temp=temp, fixcom=fixcom, eens=eens, fixatoms=fixatoms)
#~ 
   #~ def pconstraints(self):
      #~ """This removes the centre of mass contribution to the kinetic energy.
#~ 
      #~ Calculates the centre of mass momenta, then removes the mass weighted
      #~ contribution from each atom. If the ensemble defines a thermostat, then
      #~ the contribution to the conserved quantity due to this subtraction is
      #~ added to the thermostat heat energy, as it is assumed that the centre of
      #~ mass motion is due to the thermostat.
#~ 
      #~ If there is a choice of thermostats, the thermostat
      #~ connected to the centroid is chosen.
      #~ """
#~ 
      #~ if (self.fixcom):
         #~ pcom = np.zeros(3,float);
#~ 
         #~ na3 = self.beads.natoms*3
         #~ nb = self.beads.nbeads
         #~ p = depstrip(self.beads.p)
         #~ m = depstrip(self.beads.m3)[:,0:na3:3]
         #~ M = self.beads[0].M
#~ 
         #~ for i in range(3):
            #~ pcom[i] = p[:,i:na3:3].sum()
#~ 
         #~ self.eens += np.dot(pcom,pcom)/(2.0*M*nb)
#~ 
         #~ # subtracts COM _velocity_
         #~ pcom *= 1.0/(nb*M)
         #~ for i in range(3):
            #~ self.beads.p[:,i:na3:3] -= m*pcom[i]
            #~ 
      #~ if (len(self.fixatoms)>0):
         #~ for bp in self.beads.p:
            #~ m = depstrip(self.beads.m)
            #~ self.eens += 0.5*np.dot(bp[self.fixatoms*3],bp[self.fixatoms*3]/m[self.fixatoms])
            #~ self.eens += 0.5*np.dot(bp[self.fixatoms*3+1],bp[self.fixatoms*3+1]/m[self.fixatoms])
            #~ self.eens += 0.5*np.dot(bp[self.fixatoms*3+2],bp[self.fixatoms*3+2]/m[self.fixatoms])
            #~ bp[self.fixatoms*3]=0.0
            #~ bp[self.fixatoms*3+1]=0.0
            #~ bp[self.fixatoms*3+2]=0.0
#~ 
   #~ def pstep(self):
      #~ """Velocity Verlet momenta propagator."""
#~ 
      #~ self.beads.p += depstrip(self.forces.f)*(self.dt*0.5)
      #~ # also adds the bias force
      #~ self.beads.p += depstrip(self.bias.f)*(self.dt*0.5)
#~ 
   #~ def qcstep(self):
      #~ """Velocity Verlet centroid position propagator."""
#~ 
      #~ self.nm.qnm[0,:] += depstrip(self.nm.pnm)[0,:]/depstrip(self.beads.m3)[0]*self.dt
#~ 
   #~ def step(self, step=None):
      #~ """Does one simulation time step."""
#~ 
      #~ self.ptime = -time.time()
      #~ self.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.qtime = -time.time()
      #~ self.qcstep()
#~ 
      #~ self.nm.free_qstep()
      #~ self.qtime += time.time()
#~ 
      #~ self.ptime -= time.time()
      #~ self.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
#~ class NVTEnsemble(NVEEnsemble):
   #~ """Ensemble object for constant temperature simulations.
#~ 
   #~ Has the relevant conserved quantity and normal mode propagator for the
   #~ constant temperature ensemble. Contains a thermostat object containing the
   #~ algorithms to keep the temperature constant.
#~ 
   #~ Attributes:
      #~ thermostat: A thermostat object to keep the temperature constant.
#~ 
   #~ Depend objects:
      #~ econs: Conserved energy quantity. Depends on the bead kinetic and
         #~ potential energy, the spring potential energy and the heat
         #~ transferred to the thermostat.
   #~ """
#~ 
   #~ def __init__(self, dt, temp, thermostat=None, fixcom=False, eens=0.0, fixatoms=None):
      #~ """Initialises NVTEnsemble.
#~ 
      #~ Args:
         #~ dt: The simulation timestep.
         #~ temp: The system temperature.
         #~ thermostat: A thermostat object to keep the temperature constant.
            #~ Defaults to Thermostat()
         #~ fixcom: An optional boolean which decides whether the centre of mass
            #~ motion will be constrained or not. Defaults to False.
      #~ """
#~ 
      #~ super(NVTEnsemble,self).__init__(dt=dt,temp=temp, fixcom=fixcom, eens=eens, fixatoms=fixatoms)
#~ 
      #~ if thermostat is None:
         #~ self.thermostat = Thermostat()
      #~ else:
         #~ self.thermostat = thermostat
#~ 
   #~ def bind(self, beads, nm, cell, bforce, bbias, prng):
      #~ """Binds beads, cell, bforce and prng to the ensemble.
#~ 
      #~ This takes a beads object, a cell object, a forcefield object and a
      #~ random number generator object and makes them members of the ensemble.
      #~ It also then creates the objects that will hold the data needed in the
      #~ ensemble algorithms and the dependency network. Also note that the
      #~ thermostat timestep and temperature are defined relative to the system
      #~ temperature, and the the thermostat temperature is held at the
      #~ higher simulation temperature, as is appropriate.
#~ 
      #~ Args:
         #~ beads: The beads object from whcih the bead positions are taken.
         #~ nm: A normal modes object used to do the normal modes transformation.
         #~ cell: The cell object from which the system box is taken.
         #~ bforce: The forcefield object from which the force and virial are
            #~ taken.
         #~ prng: The random number generator object which controls random number
            #~ generation.
      #~ """
#~ 
      #~ super(NVTEnsemble,self).bind(beads, nm, cell, bforce, bbias, prng)
#~ 
      #~ fixdof = len(self.fixatoms)*3*self.beads.nbeads
      #~ if self.fixcom:
         #~ fixdof += 3
#~ 
#~ 
      #~ # first makes sure that the thermostat has the correct temperature, then proceed with binding it.
      #~ deppipe(self,"ntemp", self.thermostat,"temp")
      #~ deppipe(self,"dt", self.thermostat, "dt")
#~ 
      #~ #depending on the kind, the thermostat might work in the normal mode or the bead representation.
      #~ self.thermostat.bind(beads=self.beads, nm=self.nm,prng=prng,fixdof=fixdof )
#~ 
      #~ dget(self,"econs").add_dependency(dget(self.thermostat, "ethermo"))
#~ 
   #~ def step(self, step=None):
      #~ """Does one simulation time step."""
#~ 
      #~ self.ttime = -time.time()
      #~ self.thermostat.step()
      #~ self.pconstraints()
      #~ self.ttime += time.time()
#~ 
      #~ self.ptime = -time.time()
      #~ self.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.qtime = -time.time()
      #~ self.qcstep()
      #~ self.nm.free_qstep()
      #~ self.qtime += time.time()
#~ 
      #~ self.ptime -= time.time()
      #~ self.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.ttime -= time.time()
      #~ self.thermostat.step()
      #~ self.pconstraints()
      #~ self.ttime += time.time()
#~ 
   #~ def get_econs(self):
      #~ """Calculates the conserved energy quantity for constant temperature
      #~ ensemble.
      #~ """
#~ 
      #~ return NVEEnsemble.get_econs(self) + self.thermostat.ethermo
#~ 
#~ 
#~ class NPTEnsemble(NVTEnsemble):
   #~ """Ensemble object for constant pressure simulations.
#~ 
   #~ Has the relevant conserved quantity and normal mode propagator for the
   #~ constant pressure ensemble. Contains a thermostat object containing the
   #~ algorithms to keep the temperature constant, and a barostat to keep the
   #~ pressure constant.
#~ 
   #~ Attributes:
      #~ barostat: A barostat object to keep the pressure constant.
#~ 
   #~ Depend objects:
      #~ econs: Conserved energy quantity. Depends on the bead and cell kinetic
         #~ and potential energy, the spring potential energy, the heat
         #~ transferred to the beads and cell thermostat, the temperature and
         #~ the cell volume.
      #~ pext: External pressure.
   #~ """
#~ 
   #~ def __init__(self, dt, temp, pext, thermostat=None, barostat=None, fixcom=False, eens=0.0, fixatoms=None):
      #~ """Initialises NPTEnsemble.
#~ 
      #~ Args:
         #~ dt: The simulation timestep.
         #~ temp: The system temperature.
         #~ pext: The external pressure.
         #~ thermostat: A thermostat object to keep the temperature constant.
            #~ Defaults to Thermostat().
         #~ barostat: A barostat object to keep the pressure constant.
            #~ Defaults to Barostat().
         #~ fixcom: An optional boolean which decides whether the centre of mass
            #~ motion will be constrained or not. Defaults to False.
      #~ """
#~ 
      #~ super(NPTEnsemble,self).__init__(dt, temp, thermostat, fixcom=fixcom, eens=eens, fixatoms=fixatoms)
      #~ if barostat == None:
         #~ self.barostat = Barostat()
      #~ else:
         #~ self.barostat = barostat
#~ 
      #~ dset(self,"pext",depend_value(name='pext'))
      #~ if not pext is None:
         #~ self.pext = pext
      #~ else: self.pext = 0.0
#~ 
#~ 
   #~ def bind(self, beads, nm, cell, bforce, bbias, prng):
      #~ """Binds beads, cell, bforce and prng to the ensemble.
#~ 
      #~ This takes a beads object, a cell object, a forcefield object and a
      #~ random number generator object and makes them members of the ensemble.
      #~ It also then creates the objects that will hold the data needed in the
      #~ ensemble algorithms and the dependency network. Also note that the cell
      #~ thermostat timesteps and temperatures are defined relative to the system
      #~ temperature, and the the thermostat temperatures are held at the
      #~ higher simulation temperature, as is appropriate.
#~ 
      #~ Args:
         #~ beads: The beads object from whcih the bead positions are taken.
         #~ nm: A normal modes object used to do the normal modes transformation.
         #~ cell: The cell object from which the system box is taken.
         #~ bforce: The forcefield object from which the force and virial are
            #~ taken.
         #~ prng: The random number generator object which controls random number
            #~ generation.
      #~ """
#~ 
#~ 
      #~ fixdof = None
      #~ if self.fixcom:
         #~ fixdof = 3
#~ 
      #~ super(NPTEnsemble,self).bind(beads, nm, cell, bforce, bbias, prng)
      #~ self.barostat.bind(beads, nm, cell, bforce, prng=prng, fixdof=fixdof)
#~ 
#~ 
      #~ deppipe(self,"ntemp", self.barostat, "temp")
      #~ deppipe(self,"dt", self.barostat, "dt")
      #~ deppipe(self,"pext", self.barostat, "pext")
      #~ dget(self,"econs").add_dependency(dget(self.barostat, "ebaro"))
#~ 
   #~ def get_econs(self):
      #~ """Calculates the conserved energy quantity for the constant pressure
      #~ ensemble.
      #~ """
#~ 
      #~ return NVTEnsemble.get_econs(self) + self.barostat.ebaro
#~ 
   #~ def step(self, step=None):
      #~ """NPT time step.
#~ 
      #~ Note that the barostat only propagates the centroid coordinates. If this
      #~ approximation is made a centroid virial pressure and stress estimator can
      #~ be defined, so this gives the best statistical convergence. This is
      #~ allowed as the normal mode propagation is approximately unaffected
      #~ by volume fluctuations as long as the system box is much larger than
      #~ the radius of gyration of the ring polymers.
      #~ """
#~ 
      #~ self.ttime = -time.time()
      #~ self.thermostat.step()
      #~ self.barostat.thermostat.step()
      #~ self.pconstraints()
      #~ self.ttime += time.time()
#~ 
      #~ self.ptime = -time.time()
      #~ self.barostat.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.qtime = -time.time()
      #~ self.barostat.qcstep()
      #~ self.nm.free_qstep()
      #~ self.qtime += time.time()
#~ 
      #~ self.ptime -= time.time()
      #~ self.barostat.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.ttime -= time.time()
      #~ self.barostat.thermostat.step()
      #~ self.thermostat.step()
      #~ self.pconstraints()
      #~ self.ttime += time.time()
#~ 
#~ 
#~ class NSTEnsemble(NVTEnsemble):
   #~ """Ensemble object for constant pressure simulations.
#~ 
      #~ Has the relevant conserved quantity and normal mode propagator for the
      #~ constant pressure ensemble. Contains a thermostat object containing the
      #~ algorithms to keep the temperature constant, and a barostat to keep the
      #~ pressure constant.
#~ 
      #~ Attributes:
      #~ barostat: A barostat object to keep the pressure constant.
#~ 
      #~ Depend objects:
      #~ econs: Conserved energy quantity. Depends on the bead and cell kinetic
      #~ and potential energy, the spring potential energy, the heat
      #~ transferred to the beads and cell thermostat, the temperature and
      #~ the cell volume.
      #~ pext: External pressure.
      #~ """
#~ 
   #~ def __init__(self, dt, temp, stressext=None, thermostat=None, barostat=None, fixcom=False, eens=0.0, fixatoms=None):
      #~ """Initialises NSTEnsemble.
#~ 
         #~ Args:
         #~ dt: The simulation timestep.
         #~ temp: The system temperature.
         #~ stressext: The external stress.
         #~ thermostat: A thermostat object to keep the temperature constant.
         #~ Defaults to Thermostat().
         #~ barostat: A barostat object to keep the pressure constant.
         #~ Defaults to Barostat().
         #~ fixcom: An optional boolean which decides whether the centre of mass
         #~ motion will be constrained or not. Defaults to False.
         #~ """
#~ 
      #~ super(NSTEnsemble,self).__init__(dt, temp, thermostat, fixcom=fixcom, eens=eens, fixatoms=fixatoms)
      #~ if barostat == None:
         #~ self.barostat = Barostat()
      #~ else:
         #~ self.barostat = barostat
#~ 
      #~ dset(self,"stressext",depend_array(name='stressext',value=np.zeros((3,3),float)))
      #~ if not stressext is None:
         #~ self.stressext = stressext
      #~ else: self.stressext = 0.0
#~ 
#~ 
   #~ def bind(self, beads, nm, cell, bforce, bbias, prng):
      #~ """Binds beads, cell, bforce and prng to the ensemble.
#~ 
         #~ This takes a beads object, a cell object, a forcefield object and a
         #~ random number generator object and makes them members of the ensemble.
         #~ It also then creates the objects that will hold the data needed in the
         #~ ensemble algorithms and the dependency network. Also note that the cell
         #~ thermostat timesteps and temperatures are defined relative to the system
         #~ temperature, and the the thermostat temperatures are held at the
         #~ higher simulation temperature, as is appropriate.
#~ 
         #~ Args:
         #~ beads: The beads object from whcih the bead positions are taken.
         #~ nm: A normal modes object used to do the normal modes transformation.
         #~ cell: The cell object from which the system box is taken.
         #~ bforce: The forcefield object from which the force and virial are
         #~ taken.
         #~ prng: The random number generator object which controls random number
         #~ generation.
         #~ """
#~ 
#~ 
      #~ fixdof = None
      #~ if self.fixcom:
         #~ fixdof = 3
#~ 
      #~ super(NSTEnsemble,self).bind(beads, nm, cell, bforce, bbias, prng)
      #~ self.barostat.bind(beads, nm, cell, bforce, bbias, prng=prng, fixdof=fixdof)
#~ 
#~ 
      #~ deppipe(self,"ntemp", self.barostat, "temp")
      #~ deppipe(self,"dt", self.barostat, "dt")
#~ 
      #~ deppipe(self,"stressext", self.barostat, "stressext")
      #~ dget(self,"econs").add_dependency(dget(self.barostat, "ebaro"))
#~ 
   #~ def get_econs(self):
      #~ """Calculates the conserved energy quantity for the constant pressure
         #~ ensemble.
         #~ """
#~ 
      #~ return NVTEnsemble.get_econs(self) + self.barostat.ebaro
#~ 
   #~ def step(self, step=None):
      #~ """NST time step (dummy for now).
#~ 
         #~ Note that the barostat only propagates the centroid coordinates. If this
         #~ approximation is made a centroid virial pressure and stress estimator can
         #~ be defined, so this gives the best statistical convergence. This is
         #~ allowed as the normal mode propagation is approximately unaffected
         #~ by volume fluctuations as long as the system box is much larger than
         #~ the radius of gyration of the ring polymers.
         #~ """
#~ 
      #~ self.ttime = -time.time()
      #~ self.thermostat.step()
      #~ self.barostat.thermostat.step()
      #~ self.pconstraints()
      #~ self.ttime += time.time()
#~ 
      #~ self.ptime = -time.time()
      #~ self.barostat.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.qtime = -time.time()
      #~ self.barostat.qcstep()
      #~ self.nm.free_qstep()
      #~ self.qtime += time.time()
#~ 
      #~ self.ptime -= time.time()
      #~ self.barostat.pstep()
      #~ self.pconstraints()
      #~ self.ptime += time.time()
#~ 
      #~ self.ttime -= time.time()
      #~ self.barostat.thermostat.step()
      #~ self.thermostat.step()
      #~ self.pconstraints()
      #~ self.ttime += time.time()