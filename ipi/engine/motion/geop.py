"""
Contains classes for different geometry optimization algorithms.

TODO

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import depstrip, dobject
from ipi.utils.softexit import softexit
from ipi.utils.mintools import min_brent, BFGS, L_BFGS
from ipi.utils.messages import verbosity, info
from ipi.utils.counter import counter
from ipi.engine.barostats import Barostat


__all__ = ['GeopMotion']

class GeopMotion(Motion):
    """Geometry optimization class.

    Attributes:
        mode: minimization algorithm to use
        optimize_cell: boolean to decide whether cell should be optimized or not
        scale_jacobian: weighting factor for cell optimization (relatively weight importance of atomic positions and cell)
        biggest_step: max allowed step size for BFGS/L-BFGS
        old_force: force on previous step
        old_direction_cgsd: move direction on previous step in CG/SD
        invhessian_bfgs: stored inverse Hessian matrix for BFGS
        ls_options:
        {tolerance: energy tolerance for exiting minimization algorithm
        iter: maximum number of allowed iterations for minimization algorithm for each MD step
        step: initial step size for steepest descent and conjugate gradient
        adaptive: T/F adaptive step size for steepest descent and conjugate
                gradient}
        tolerances:
        {energy: change in energy tolerance for ending minimization
        force: force/change in force tolerance foe ending minimization
        position: change in position tolerance for ending minimization}
        corrections_lbfgs: number of corrections to be stored for L-BFGS
        qlist_lbfgs: list of previous positions (x_n+1 - x_n) for L-BFGS. Number of entries = corrections_lbfgs
        glist_lbfgs: list of previous gradients (g_n+1 - g_n) for L-BFGS. Number of entries = corrections_lbfgs
    """

    def __init__(self, fixcom=False, fixatoms=None,
                 mode="lbfgs",
                 optimize_cell=True,
                 scale_jacobian=1.0,
                 biggest_step=100.0,
                 old_force=np.zeros(0, float),
                 old_direction_cgsd=np.zeros(0, float),
                 invhessian_bfgs=np.eye(0, 0, 0, float),
                 ls_options={"tolerance": 1e-6, "iter": 100, "step": 1e-3, "adaptive": 1.0},
                 tolerances={"energy": 1e-8, "force": 1e-8, "position": 1e-8},
                 corrections_lbfgs=5,
                 qlist_lbfgs=np.zeros(0, float),
                 glist_lbfgs=np.zeros(0, float)):
        """Initialises GeopMotion.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
        """

        super(GeopMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)
        
        # Optimization Options

        self.mode = mode
        self.optimize_cell = optimize_cell
        self.scale_jacobian = scale_jacobian
        self.big_step = biggest_step
        self.old_f = old_force
        self.old_d = old_direction_cgsd
        self.invhessian = invhessian_bfgs
        self.ls_options = ls_options
        self.tolerances = tolerances
        self.corrections = corrections_lbfgs
        self.qlist = qlist_lbfgs
        self.glist = glist_lbfgs
        self.qcell = np.zeros(0, dtype=float)
        
        # Classes for minimization routines
        self.optype = mode
        if self.optype == "bfgs":
            self.optimizer = BFGSOptimizer()
        elif self.optype == "lbfgs":
            self.optimizer = LBFGSOptimizer()
        elif self.optype == "sd":
            self.optimizer = SDOptimizer()
        elif self.optype == "cg":
            self.optimizer = CGOptimizer()
        else:
            self.optimizer = DummyOptimizer()
        
        
    def bind(self, ens, beads, nm, cell, bforce, prng):
        """Binds beads, cell, bforce and prng to GeopMotion
        
            Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            prng: The random number generator object which controls random number generation.
        """

        super(GeopMotion,self).bind(ens, beads, nm, cell, bforce, prng)
        # Binds optimizer (call bind function from DummyOptimizer)
        self.optimizer.bind(self)
       
        
    def step(self, step=None):
        self.optimizer.step(step)

class Mapper(object):

    """Parent class for function evaluations during optimization.
       Attributes:
       x0: initial position
       d: move direction"""

    def __init__(self):
        self.x0 = None
        self.d = None
        self.strain = None
        self.oldcell = None
   
    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)
        self.dqcell = dumop.qcell.copy()
        self.optimize_cell = dumop.optimize_cell
        self.jacobian = dumop.jacobian
        self.norm_stress = dumop.norm_stress
    
    def transform(self, x):

        """Function to transform the stress and the strain into the units of forces and atomic positions; also updates cell."""
        natoms = self.dbeads.natoms
        self.dbeads.q[-1,:] = x[0:natoms*3]      # atomic positions
        self.strain = x[natoms*3:]/self.jacobian # strain vector
        # we use the Voight strain tensor.
        #Lower part is set to zero because i-Pi expects lower triangular cell matrix. 
        #Should give the same result (tested with ase) but takes much more optimization steps.
        epsilon =  np.array([[1.0 + self.strain[0], 0.5 * self.strain[5], 0.5 * self.strain[4]],
                        [0.0 ,       1.0 + self.strain[1],         0.5 * self.strain[3]],
                        [0.0,         0.0 ,                   1.0 + self.strain[2]]])
        self.dcell.h = np.dot(self.oldcell, epsilon)     #update cell
        stress = self.dforces.vir.flatten()/self.dcell.V #stress tensor
        stress_force = np.array([stress[0], stress[4], stress[8],
                                 stress[5], stress[2], stress[1]])*self.norm_stress #scaled stress vector
        return stress_force
        
        
class LineMapper(Mapper):
    
    """Creation of the one-dimensional function that will be minimized.
    Used in steepest descent and conjugate gradient minimizers.

    Attributes:
        x0: initial position
        d: move direction
    """

    def __init__(self):
        super(LineMapper, self).__init__()

    def bind(self, dumop):
        super(LineMapper, self).bind(dumop)

    def set_dir(self, x0, mdir):
        self.x0 = x0.copy()
        self.d =  mdir.copy() / np.sqrt(np.dot(mdir.flatten(), mdir.flatten()))
        if self.x0.shape != self.d.shape:
            raise ValueError("Incompatible shape of initial value and displacement direction")

    def __call__(self, x):
        """ computes energy and gradient for optimization step
            determines new position (x0+d*x)"""
        
        if self.optimize_cell == True:
            self.dqcell = self.x0 + self.d * x      # new positions
            stress_force = self.transform(self.dqcell)
            natoms = self.dbeads.natoms
            forces = np.zeros((natoms + 2)*3, float)
            forces[0:natoms*3] = depstrip(self.dforces.f.flatten())
            forces[natoms*3:]  = stress_force
        else:
            self.dbeads.q = self.x0 + self.d * x    # new positions
            forces = depstrip(self.dforces.f.flatten())
        
        e = self.dforces.pot                     # Energy
        g = - np.dot(forces, self.d.flatten())   # Gradient
        counter.count()                          # counts number of function evaluations
        return e, g

class GradientMapper(Mapper):
       
    """Creation of the multi-dimensional function that will be minimized.
    Used in the BFGS and L-BFGS minimizers.

    Attributes:
        x0: initial position
        d: move direction
        xold: previous position
    """
    
    def __init__(self):
        super(GradientMapper, self).__init__()
        self.xold = None
        
    def bind(self, dumop):
         super(GradientMapper, self).bind(dumop)
        
    def __call__(self,x):
        """computes energy and gradient for optimization step"""
        
        if self.optimize_cell == True:
            stress_force = self.transform(x)
            natoms = self.dbeads.natoms
            forces = np.zeros((natoms + 2)*3, float)
            forces[0:natoms*3] = depstrip(self.dforces.f.flatten())
            forces[natoms*3:]  = stress_force
        else:
            self.dbeads.q = x  # new positions
            forces = self.dforces.f
        e = self.dforces.pot   # Energy
        g = - forces           # Gradient
        counter.count()        # counts number of function evaluations
        return e, g

class DummyOptimizer(dobject):
    """ Dummy class for all optimization classes """
    
    def __init__(self):
        """initialises object for LineMapper (1-d function) and for GradientMapper (multi-dimensional function) """
        self.m = Mapper()
        self.lm = LineMapper()
        self.gm = GradientMapper()
        
    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass
        
    def bind(self, geop):
        """ 
        bind optimization options and call bind function of LineMapper and GradientMapper (get beads, cell,forces)
        check whether force size, direction size and inverse Hessian size from previous step match system size
        """
        
        self.ls_options = geop.ls_options   
        self.tolerances = geop.tolerances
        self.mode = geop.mode               
        self.big_step = geop.big_step       
        self.old_f = geop.old_f
        self.old_d = geop.old_d
        self.invhessian = geop.invhessian   
        self.corrections = geop.corrections 
        self.qlist = geop.qlist             
        self.glist = geop.glist             
        self.beads = geop.beads
        self.cell = geop.cell
        self.forces = geop.forces
        self.fixcom = geop.fixcom
        self.fixatoms = geop.fixatoms
        self.qcell = geop.qcell
        self.optimize_cell = geop.optimize_cell
        self.scale_jacobian = geop.scale_jacobian

        #Jacobians to scale the stress (norm_stress) and the strain (jacobian) with the system size and to give them the same units as atomic positions/forces
        # method: http://scitation.aip.org/content/aip/journal/jcp/136/7/10.1063/1.3684549

        self.jacobian = self.scale_jacobian * self.cell.V**(1.0/3.0)*self.beads.natoms**(1.0/6.0)
        self.norm_stress = -self.cell.V/self.jacobian

        self.lm.bind(self)
        self.gm.bind(self)

        if self.optimize_cell == True:
            #vector that contains atomic positions and lattice parameters
            self.qcell = np.zeros(3*(self.beads.natoms+2), float)
        else:
            #vector that contains only atomic positions
            self.qcell = np.zeros(3*(self.beads.natoms), float)
        
        if self.old_f.shape != self.qcell.size:
            if self.old_f.size == 0:
                self.old_f = np.zeros(self.qcell.size, float)
            else:
                raise ValueError("Conjugate gradient force size does not match system size")
        if self.old_d.size != self.qcell.size:
            if self.old_d.size == 0:
                self.old_d = np.zeros(self.qcell.size, float)
            else:
                raise ValueError("Conjugate gradient direction size does not match system size")
        if self.invhessian.size != ((self.qcell.size) * (self.qcell.size)):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(self.qcell.size, self.qcell.size, 0, float)
            else:
                raise ValueError("Inverse Hessian size does not match system size")
        
    def exitstep(self, fx, u0, x, forces):
        """ Exits the simulation step. Computes time, checks for convergence. """
        
        info(" @GEOP: Updating bead positions", verbosity.debug)
        
        self.qtime += time.time()

        # Determine conditions for converged relaxation
        if (np.absolute(fx - u0) / (self.qcell.size/3) <= self.tolerances["energy"])\
                and ((np.amax(np.absolute(forces)) <= self.tolerances["force"])
                   # or (np.sqrt(np.dot(forces.flatten() - self.old_f.flatten(),
                       # forces.flatten() - self.old_f.flatten())) == 0.0))\
            #changed by Sophie
          or (np.sqrt(np.dot(forces.flatten() - self.old_f.flatten(), forces.flatten() - self.old_f.flatten())) <= 1e-10))\
                and (x <= self.tolerances["position"]):
            info("Total number of function evaluations: %d" % counter.func_eval, verbosity.debug)
            softexit.trigger("Geometry optimization converged. Exiting simulation")
            
    def update_cell_q(self):
        """Function to update the positions and the cell in a cell optimization."""

        # number of atoms
        natoms = self.beads.natoms
        self.beads.q[-1,:] = self.qcell[0:natoms*3]
        self.m.strain = self.qcell[natoms*3:]/self.jacobian #strain vector
        # we use the Voight strain tensor.
        #Lower part is set to zero because i-Pi expects lower triangular cell matrix. 
        #Should give the same result (tested with ase) but takes much more optimization steps.
        epsilon =  np.array([[1.0 + self.m.strain[0], 0.5 * self.m.strain[5], 0.5 * self.m.strain[4]],
                        [0.0 ,       1.0 + self.m.strain[1],         0.5 * self.m.strain[3]],
                        [0.0,         0.0 ,                   1.0 + self.m.strain[2]]])
        return epsilon
       
    def transform_force(self):
        """Function to transform stress into the units of forces and put scaled stress and forces in one vector. """

        # number of atoms
        natoms = self.beads.natoms
        # new force vector with forces and transformed stress tensor
        forces = np.zeros((natoms+2)*3, float)
        forces[0:natoms*3] = self.forces.f.flatten()
        stress = self.forces.vir.flatten()/self.cell.V
        stress_force = np.array([stress[0], stress[4], stress[8],
                  stress[5], stress[2], stress[1]])*self.norm_stress
        forces[natoms*3:] = stress_force
        return forces
    
class BFGSOptimizer(DummyOptimizer):
    """ BFGS Minimization """

    def step(self, step=None):
        """ Does one simulation time step.
            Attributes:
            ptime: The time taken in updating the velocities.
            qtime: The time taken in updating the positions.
            ttime: The time taken in applying the thermostat steps.
        """

        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        info("\nMD STEP %d" % step, verbosity.debug)

        # Initialize approximate Hessian inverse to the identity and direction
        # to the steepest descent direction

        natoms = self.beads.natoms
        if step == 0:   # or np.sqrt(np.dot(self.gm.d, self.gm.d)) == 0.0: this part for restarting at claimed minimum (optional)
            info(" @GEOP: Initializing BFGS", verbosity.debug)
            if self.optimize_cell == True:
                # store actual position to previous position, initial strain is zero
                self.gm.xold = np.zeros((natoms+2)*3, float)
                self.gm.xold[0:natoms*3] = self.beads.q.copy()
                self.m.strain = np.zeros(6, float)
                # store cell
                self.gm.oldcell = self.cell.h.copy()
                #force vector containing atomic forces and scaled stress
                forces = self.transform_force()
                # gradient direction
                self.gm.d = depstrip(forces) / np.sqrt(np.dot(forces, forces))
            else:
                # store actual position to previous position
                self.gm.xold = self.beads.q.copy()
                # gradient direction
                self.gm.d = depstrip(self.forces.f.flatten()) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
        
        if self.optimize_cell == True:
            #new vector containing atomic positions and transformed strain, initial strain is zero
            self.qcell[0:natoms*3] = self.beads.q.copy()
            self.qcell[natoms*3:] = 0.0
            #force vector containing atomic forces and scaled stress
            forces = self.transform_force()
        else:
            self.qcell = self.beads.q.copy()
            forces = self.forces.f.flatten()

        # Current energy and forces
        u0 = self.forces.pot.copy()
        du0 = - forces

        # Store previous forces
        self.old_f[:] = forces

        # Do one iteration of BFGS, return new point, function value,
        # move direction, and current Hessian to use for next iteration
        self.qcell, fx, self.gm.d, self.invhessian = BFGS(self.qcell,
                self.gm.d, self.gm, fdf0=(u0, du0), invhessian=self.invhessian,
                big_step=self.big_step, tol=self.ls_options["tolerance"],
                                                            itmax=self.ls_options["iter"])

        # x = current position - previous position; use for exit tolerance
        x = np.amax(np.absolute(np.subtract(self.qcell, self.gm.xold)))
        
        if self.optimize_cell == True:
            # update cell, atomic positions
            epsilon = self.update_cell_q()
            self.cell.h = np.dot(self.gm.oldcell, epsilon)
            self.gm.oldcell = self.cell.h.copy()
            # this is to store the old strain as zero (for gm.xold)
            self.qcell[natoms*3:] = 0.0
            forces = self.transform_force()
        else:
            self.beads.q = self.qcell
            forces = self.forces.f.flatten()
        
        # Store old position
        self.gm.xold[:] = self.qcell
        
        # Exit simulation step
        self.exitstep(fx, u0, x, forces)

class LBFGSOptimizer(DummyOptimizer):
    """ L-BFGS Minimization """
    
    def step(self, step=None):
        """ Does one simulation time step 
            Attributes:
            ptime: The time taken in updating the velocities.
            qtime: The time taken in updating the positions.
            ttime: The time taken in applying the thermostat steps.
        """

        natoms = self.beads.natoms

        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        info("\nMD STEP %d" % step, verbosity.debug)
        
        # Initialize approximate Hessian inverse to the identity and direction
        # to the steepest descent direction
        if step == 0:   # or np.sqrt(np.dot(self.gm.d, self.gm.d)) == 0.0: <-- this part for restarting at claimed minimum (optional)
            info(" @GEOP: Initializing L-BFGS", verbosity.debug)
            if self.optimize_cell == True:
                # store actual position to previous position, initial strain is zero
                self.gm.xold = np.zeros((natoms+2)*3, float)
                self.gm.xold[0:natoms*3] = self.beads.q.copy()
                self.m.strain = np.zeros(6, float)
                # store cell
                self.gm.oldcell = self.cell.h.copy()
                # force vector containing atomic forces and scaled stress
                forces = self.transform_force()
            else:
                # store actual position to previous position
                self.gm.xold = self.beads.q.copy()
                forces = self.forces.f.flatten()
            
            # gradient direction
            self.gm.d = depstrip(forces) / np.sqrt(np.dot(forces.flatten(), forces.flatten()))
            # Initialize lists of previous positions and gradient
            self.qlist = np.zeros((self.corrections, len(self.qcell)))
            self.glist = np.zeros((self.corrections, len(self.qcell)))

        if self.optimize_cell == True:
            #new vector containing atomic positions and transformed strain, initial strain is zero
            self.qcell[0:natoms*3] = self.beads.q.copy()
            self.qcell[natoms*3:] = 0.0
            # force vector containing atomic forces and scaled stress
            forces = self.transform_force()
        else:
            self.qcell = self.beads.q.copy()
            forces = self.forces.f.flatten()

        # Current energy and force
        u0, du0 = (self.forces.pot.copy(), - forces)

        # Store previous forces
        self.old_f[:] = forces

        # Do one iteration of L-BFGS, return new point, function value,
        # move direction, and current Hessian to use for next iteration
        self.qcell, fx, self.gm.d, self.qlist, self.glist = L_BFGS(self.qcell,
                self.gm.d, self.gm, self.qlist, self.glist,
                fdf0=(u0, du0), big_step=self.big_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"],
                m=self.corrections, k=step)

        info(" @GEOP: Updated position list", verbosity.debug)
        info(" @GEOP: Updated gradient list", verbosity.debug)

        if self.optimize_cell == True:
            # update positions, cell
            epsilon = self.update_cell_q()
            self.cell.h = np.dot(self.gm.oldcell, epsilon)
            self.gm.oldcell = self.cell.h.copy()
            # this is to store the old strain as zero (for gm.xold)
            self.qcell[natoms*3:]= 0.0
            forces = self.transform_force()
        else:
            # update positions
            self.beads.q = self.qcell
            forces = self.forces.f

        # x = current position - old position. Used for convergence tolerance
        x = np.amax(np.absolute(np.subtract(self.qcell, self.gm.xold)))
        
        # Store old position
        self.gm.xold[:] = self.qcell
        
        # Exit simulation step
        self.exitstep(fx, u0, x, forces)
         
class SDOptimizer(DummyOptimizer):
    """
    Steepest descent minimization
    gradf1 = force at current atom position
    dq1 = direction of steepest descent
    dq1_unit = unit vector of dq1
    """

    def step(self, step=None):
        """ Does one simulation time step 
            Attributes:
            ptime: The time taken in updating the velocities.
            qtime: The time taken in updating the positions.
            ttime: The time taken in applying the thermostat steps.
        """
        
        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        info("\nMD STEP %d" % step, verbosity.debug)
        
        natoms = self.beads.natoms
        if step == 0:
            # store cell
            self.lm.oldcell = self.cell.h.copy()

        if self.optimize_cell == True:
            # vector containing atomic positions and scaled strain, initial strain is zero
            self.qcell[0:natoms*3] = self.beads.q.copy()
            self.qcell[natoms*3:] = 0.0
            # vector containing atomic forces and scaled stress
            forces = self.transform_force()
        else:
            self.qcell = self.beads.q.copy()
            forces = self.forces.f

        gradf1 = dq1 = depstrip(forces)

        # Move direction for steepest descent
        dq1_unit = dq1 / np.sqrt(np.dot(gradf1.flatten(), gradf1.flatten()))
        info(" @GEOP: Determined SD direction", verbosity.debug)
       
        # Store force and direction for next CG step
        self.old_d[:] = dq1
        self.old_f[:] = gradf1
        
        if len(self.fixatoms) > 0:
            for dqb in dq1_unit:
                dqb[self.fixatoms*3] = 0.0
                dqb[self.fixatoms*3+1] = 0.0
                dqb[self.fixatoms*3+2] = 0.0
        
        self.lm.set_dir(depstrip(self.qcell), dq1_unit)
        
        # Reuse initial value since we have energy and forces already
        u0, du0 = (self.forces.pot.copy(), np.dot(depstrip(forces.flatten()), dq1_unit.flatten()))

        # Do one SD iteration; return positions and energy
        (x, fx) = min_brent(self.lm, fdf0=(u0, du0), x0=0.0,
                    tol=self.ls_options["tolerance"],
                    itmax=self.ls_options["iter"], init_step=self.ls_options["step"])
        # Automatically adapt the search step for the next iteration.
        # Relaxes better with very small step --> multiply by factor of 0.1 or 0.01
        self.ls_options["step"] = 0.1 * x * self.ls_options["adaptive"] + (1 - self.ls_options["adaptive"]) * self.ls_options["step"]

        if self.optimize_cell == True:
            # update positions, cell
            self.qcell += dq1_unit * x
            epsilon = self.update_cell_q()
            self.cell.h = np.dot(self.lm.oldcell, epsilon)
            self.lm.oldcell = self.cell.h.copy()
            forces = self.transform_force()
        else:
            # update positions
            self.beads.q += dq1_unit * x
            forces = self.forces.f
        
        # Exit simulation step
        self.exitstep(fx, u0, x, forces)
        

class CGOptimizer(DummyOptimizer):
    """
    Conjugate gradient, Polak-Ribiere
    gradf1: force at current atom position
    gradf0: force at previous atom position
    dq1 = direction to move
    dq0 = previous direction
    dq1_unit = unit vector of dq1
    """

    def step(self, step=None):
        """Does one simulation time step 
           Attributes:
           ptime: The time taken in updating the velocities.
           qtime: The time taken in updating the positions.
           ttime: The time taken in applying the thermostat steps.
        """
        
        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        info("\nMD STEP %d" % step, verbosity.debug)
        natoms = self.beads.natoms

        if self.optimize_cell == True:
            # vector containing atomic positions and scaled strain, initial strain is zero
            self.qcell[0:natoms*3] = self.beads.q.copy()
            self.qcell[natoms*3:] = 0.0
            # vector containing atomic forces and scaled stress
            forces = self.transform_force()
        else:
            self.qcell = self.beads.q.copy()
            forces = self.forces.f

        if step == 0:
            gradf1 = dq1 = depstrip(forces)
            # store cell
            self.lm.oldcell = self.cell.h.copy()

            # Move direction for 1st conjugate gradient step
            dq1_unit = dq1 / np.sqrt(np.dot(gradf1.flatten(), gradf1.flatten()))
            info(" @GEOP: Determined SD direction", verbosity.debug)
    
        else:
        
            gradf0 = self.old_f
            dq0 = self.old_d
            gradf1 = depstrip(forces)
            beta = np.dot((gradf1.flatten() - gradf0.flatten()), gradf1.flatten()) / (np.dot(gradf0.flatten(), gradf0.flatten()))
            dq1 = gradf1 + max(0.0, beta) * dq0
            dq1_unit = dq1 / np.sqrt(np.dot(dq1.flatten(), dq1.flatten()))
            info(" @GEOP: Determined CG direction", verbosity.debug)

        # Store force and direction for next CG step
        self.old_d[:] = dq1
        self.old_f[:] = gradf1

        if len(self.fixatoms) > 0:
            for dqb in dq1_unit:
                dqb[self.fixatoms*3] = 0.0
                dqb[self.fixatoms*3+1] = 0.0
                dqb[self.fixatoms*3+2] = 0.0

        self.lm.set_dir(depstrip(self.qcell), dq1_unit)

        # Reuse initial value since we have energy and forces already
        u0, du0 = (self.forces.pot.copy(), np.dot(depstrip(forces), dq1_unit.flatten()))

        # Do one CG iteration; return positions and energy
        (x, fx) = min_brent(self.lm, fdf0=(u0, du0), x0=0.0,
                    tol=self.ls_options["tolerance"],
                    itmax=self.ls_options["iter"], init_step=self.ls_options["step"])

        # Automatically adapt the search step for the next iteration.
        # Relaxes better with very small step --> multiply by factor of 0.1 or 0.01
        self.ls_options["step"] = 0.1 * x * self.ls_options["adaptive"] + (1 -     self.ls_options["adaptive"]) * self.ls_options["step"]
        
        if self.optimize_cell == True:
            #update positions, cell
            self.qcell += dq1_unit * x
            epsilon = self.update_cell_q()
            self.cell.h = np.dot(self.lm.oldcell, epsilon)
            self.lm.oldcell = self.cell.h.copy()
            forces = self.transform_force()
        else:
            #update positions
            self.beads.q += dq1_unit * x
            forces = self.forces.f
        
        # Exit simulation step
        self.exitstep(fx, u0, x, forces)
