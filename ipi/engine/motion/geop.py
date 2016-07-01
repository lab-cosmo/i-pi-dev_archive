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
            self.optimizer = BFGSCellOptimizer()
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
        
class LineMapper(object):
    
    """Creation of the one-dimensional function that will be minimized.
    Used in steepest descent and conjugate gradient minimizers.

    Attributes:
        x0: initial position
        d: move direction
    """
    
    def __init(self):
        self.x0 = self.d = None
        
    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)
    
    def set_dir(self, x0, mdir):
        self.x0 = x0.copy()
        self.d = mdir.copy() / np.sqrt(np.dot(mdir.flatten(), mdir.flatten()))
        if self.x0.shape != self.d.shape:
            raise ValueError("Incompatible shape of initial value and displacement direction")
            
    def __call__(self, x):
        """ computes energy and gradient for optimization step
            determines new position (x0+d*x)"""
        
        self.dbeads.q = self.x0 + self.d * x
        e = self.dforces.pot   # Energy
        g = - np.dot(depstrip(self.dforces.f).flatten(), self.d.flatten())   # Gradient
        counter.count()      # counts number of function evaluations
        return e, g
    
class GradientMapper(object):
       
    """Creation of the multi-dimensional function that will be minimized.
    Used in the BFGS and L-BFGS minimizers.

    Attributes:
        x0: initial position
        d: move direction
        xold: previous position
    """
    
    def __init__(self):
        self.x0 = None
        self.d = None
        self.xold = None
        
    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)
        
    def __call__(self,x):
        """computes energy and gradient for optimization step"""
        
        self.dbeads.q = x
        e = self.dforces.pot   # Energy
        g = - self.dforces.f   # Gradient
        counter.count()        # counts number of function evaluations
        return e, g

class GradientCellMapper(object):
       
    """Creation of the multi-dimensional function that will be minimized.
    Used in the BFGS and L-BFGS minimizers. Testclass for minimization that also includes the cell parameters.

    Attributes:
        x0: initial position
        d: move direction
        xold: previous position
    """
    
    def __init__(self):
        self.x0 = None
        self.d = None
        self.xold = None
        self.strain = None
        self.oldcell = None
        
    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)
        self.dqcell = dumop.qcell.copy()
        self.jacobian = dumop.jacobian
        self.norm_stress = dumop.norm_stress
    
    def transform(self, x):
        '''gets a vector containing the new atomic positions and the transformed strain. Transforms it back to strain and updates positions and cell parameters'''
        natoms = self.dbeads.natoms
        self.dqcell = x
        self.dbeads.q[-1,:] = x[0:natoms*3]
        eps_vec = x[natoms*3:]/self.jacobian
        #eps = np.reshape(eps_vec, (3,3))
        #print eps
        eps = np.array([[1.0 + eps_vec[0], 0.5 * eps_vec[5], 0.5 * eps_vec[4]],
                        [0.5 * eps_vec[5], 1.0 + eps_vec[1], 0.5 * eps_vec[3]],
                        [0.5 * eps_vec[4], 0.5 * eps_vec[3], 1.0 + eps_vec[2]]])
        #unit = np.eye(3,dtype=float)
        #self.dcell.h = np.dot(self.oldcell, unit + eps)
        self.dcell.h = np.dot(self.oldcell, eps)
        M = np.linalg.solve(self.oldcell, self.dcell.h)
        #print 'M', M
        #qresh=np.reshape(x[0:natoms*3],(natoms,3))
        #for i in range(natoms):
        #    qresh[i,:]=np.dot(qresh[i,:],M)
        #self.dbeads.q[-1,:] = qresh.flatten()
        #self.dqcell=x
        #self.dqcell[0:natoms*3] = self.dbeads.q[-1,:]
        self.strain = eps_vec
        
    def transformcell(self, x):
        self.oldcell = self.dcell.h.copy()
        natoms = self.dbeads.natoms
        self.dqcell = x
        eps_vec = x/self.jacobian
        eps = np.reshape(eps_vec, (3,3))
        unit = np.eye(3,dtype=float)
        self.dcell.h = np.dot(self.oldcell, unit + eps)
        #self.dcell.h = np.dot(self.oldcell, eps)
        self.strain[:] = eps_vec
        print 'Strain during BFGS steps', self.strain
            
    def __call__(self,x):
        """computes energy and gradient for optimization step and creates vector containing atomic forces and transformed stress """
        
        self.transform(x)
        natoms = self.dbeads.natoms
        e = self.dforces.pot   # Energy
        #print 'energy', e
        g = np.zeros((natoms + 2)*3, float)
        g[0:natoms*3] = - self.dforces.f
        stress = self.dforces.vir.flatten()
        stressmat = np.array([stress[0], stress[4], stress[8],
                        stress[5], stress[2], stress[1]])
        g[natoms*3:] = -stressmat*self.norm_stress/self.dcell.V
        counter.count()        # counts number of function evaluations
        return e, g          

class DummyOptimizer(dobject):
    """ Dummy class for all optimization classes """
    
    def __init__(self):
        """initialises object for LineMapper (1-d function) and for GradientMapper (multi-dimensional function) """
        
        self.lm = LineMapper()
        self.gm = GradientMapper()
        self.gmc = GradientCellMapper()
        
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
        #Jacobians to scale the stress (norm_stress) and the strain (jacobian) with the system size and giving them the same units as atomic positions/forces
        self.jacobian = self.cell.V**(1.0/3.0)*self.beads.natoms**(1.0/6.0)
        self.norm_stress = -self.cell.V/self.jacobian


        self.lm.bind(self)
        self.gm.bind(self)
        #vector that contains atomic positions and cell parameters
        self.gmc.bind(self)
        self.qcell = np.zeros(3*(self.beads.natoms+2), float)
        #self.qcell=np.zeros(6,float)
        '''
        if self.old_f.shape != self.beads.q.size:
            if self.old_f.size == 0:
                self.old_f = np.zeros(self.beads.q.size, float)
            else:
                raise ValueError("Conjugate gradient force size does not match system size")
        if self.old_d.size != self.beads.q.size:
            if self.old_d.size == 0:
                self.old_d = np.zeros(self.beads.q.size, float)
            else:
                raise ValueError("Conjugate gradient direction size does not match system size")
        if self.invhessian.size != (self.beads.q.size * self.beads.q.size):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(self.beads.q.size, self.beads.q.size, 0, float)
            else:
                raise ValueError("Inverse Hessian size does not match system size")
        ''' 
        
        if self.old_f.shape != self.beads.q.size+6:
            if self.old_f.size == 0:
                self.old_f = np.zeros(self.beads.q.size+6, float)
            else:
                raise ValueError("Conjugate gradient force size does not match system size")
        if self.old_d.size != self.beads.q.size+6:
            if self.old_d.size == 0:
                self.old_d = np.zeros(self.beads.q.size+6, float)
            else:
                raise ValueError("Conjugate gradient direction size does not match system size")
        if self.invhessian.size != ((self.beads.q.size+6) * (self.beads.q.size+6)):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(self.beads.q.size+6, self.beads.q.size+6, 0, float)
            else:
                raise ValueError("Inverse Hessian size does not match system size")
        '''
        if self.old_f.shape != 9:
            if self.old_f.size == 0:
                self.old_f = np.zeros(9, float)
            else:
                raise ValueError("Conjugate gradient force size does not match system size")
        if self.old_d.size != 9:
            if self.old_d.size == 0:
                self.old_d = np.zeros(9, float)
            else:
                raise ValueError("Conjugate gradient direction size does not match system size")
        if self.invhessian.size != (9 * 9):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(9, 9, 0, float)
            else:
                raise ValueError("Inverse Hessian size does not match system size")
        '''
    def exitstep(self, fx, u0, x):
        """ Exits the simulation step. Computes time, checks for convergence. """
        
        info(" @GEOP: Updating bead positions", verbosity.debug)
        
        self.qtime += time.time()
        
        # Determine conditions for converged relaxation
        if ((fx - u0) / self.beads.natoms <= self.tolerances["energy"])\
                and ((np.amax(np.absolute(self.forces.f)) <= self.tolerances["force"])
                    or (np.sqrt(np.dot(self.forces.f.flatten() - self.old_f.flatten(),
                        self.forces.f.flatten() - self.old_f.flatten())) == 0.0))\
                and (x <= self.tolerances["position"]):
            info("Total number of function evaluations: %d" % counter.func_eval, verbosity.debug)
            softexit.trigger("Geometry optimization converged. Exiting simulation")
            
    def exitstep2(self, fx, u0, x, forces):
        """ Exits the simulation step. Computes time, checks for convergence. """
        
        info(" @GEOP: Updating bead positions", verbosity.debug)
        
        self.qtime += time.time()
        # Determine conditions for converged relaxation
        print '1', np.amax(np.absolute(forces))
        print '2', (np.sqrt(np.dot(forces - self.old_f.flatten(),
                        forces - self.old_f.flatten())))
        if ((fx - u0) / (self.beads.natoms+2) <= self.tolerances["energy"])\
                and (x <= self.tolerances["position"]):
            print self.beads.q*0.52917721, self.cell.h*0.52917721
            info("Total number of function evaluations: %d" % counter.func_eval, verbosity.debug)
            softexit.trigger("Geometry optimization converged. Exiting simulation")
            '''                and ((np.amax(np.absolute(forces)) <= self.tolerances["force"])
                    or (np.sqrt(np.dot(forces - self.old_f.flatten(),
                        forces - self.old_f.flatten())) == 0.0))\
                        '''
            
    def exitstepsimple(self, fx, u0, x, forces):
        """ Exits the simulation step. Computes time, checks for convergence. """
        
        info(" @GEOP: Updating bead positions", verbosity.debug)
        
        self.qtime += time.time()
        # Determine conditions for converged relaxation
        print 'TOLERANCE', np.amax(np.absolute(forces))
        #print forces
        if  ((np.amax(np.absolute(forces)) <= self.tolerances["force"])):
            info("Total number of function evaluations: %d" % counter.func_eval, verbosity.low)
            softexit.trigger("Geometry optimization converged. Exiting simulation")
            


'''
class TestVirial(DummyOptimizer):
    
    def step(self, step=None):
        u = self.forces.pot.copy()
        v = self.forces.vir
        print 'virial', self.forces.vir/self.cell.V
        print self.cell.h
        for i in range(0,3):
            for j in range(0,3):
                hdel = self.cell.h.copy()
                hdel[i,j]=hdel[i,j]+0.0001
                udel, forces = self.gm(self.beads.q, hdel)
                hminus = self.cell.h.copy()
                hminus[i,j]=hminus[i,j]-0.0001
                uminus, forces = self.gm(self.beads.q, hminus)
                vtest = (udel-uminus)/0.0002
                
                print 'vir',i,j, 	np.dot(self.beads.q[-1][i:3:3*self.beads.natoms],self.forces.f[-1][j:3:3*self.beads.natoms])/self.cell.V
                
        softexit.trigger('Blablabla')
       



class TestBFGS(DummyOptimizer):
    
    def step(self, step=None):
        
        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()
        
        natoms = self.beads.natoms
        print self.forces.pot
        print 'step', step
        if step == 0:
            #forces = np.zeros((natoms+2)*3, float)
            #forces[natoms*3] = 1.0*self.norm_stress
            forces = np.zeros((natoms+2)*3, float)
            forces[0:natoms*3]=self.forces.f
            stress = self.forces.vir.flatten()/self.cell.V
            stressnew = np.array([[stress[0], stress[4], stress[8]],
                        [stress[5], stress[2], stress[1]]])
            forces[natoms*3:] = stressnew.flatten()*self.norm_stress
            g = -forces/ np.sqrt(np.dot(forces, forces))
            #forces = np.zeros(6, float)
            #forces[0] = 1.0*self.norm_stress
            #forces[1] = -1.0*self.norm_stress
            #forces[5] = 2.0*self.norm_stress
            #g = -forces/np.sqrt(np.dot(forces, forces))
            self.gmc.strain = np.zeros(6, float)
            #self.gmc.xold = np.zeros((natoms+3)*3, float)
            #self.gmc.xold[0:natoms*3] = self.beads.q.copy()
            self.gmc.xold = np.zeros((natoms+2)*3, float)
            self.gmc.oldcell = self.cell.h.copy()
            
        print'strain', self.gmc.strain
        
        if step > 0:
            forces = np.zeros((natoms+2)*3, float)
            forces[0:natoms*3]=self.forces.f
            #forces[natoms*3:] = self.forces.vir.copy()/self.cell.V
            #forces = np.zeros(6, float)
            #forces[0] = 1.0*self.norm_stress
            #forces[4] = -1.0*self.norm_stress
            #forces[1] = 2.0*self.norm_stress
            stress = self.forces.vir.flatten()/self.cell.V
            stressnew = np.array([[stress[0], stress[4], stress[8]],
                        [stress[5], stress[2], stress[1]]])
            forces[natoms*3:] = stressnew.flatten()*self.norm_stress
            #forces = stress.flatten()*self.norm_stress
            g = -forces/ np.sqrt(np.dot(forces, forces))
            print 'dotforces', np.sqrt(np.dot(forces, forces))
            if step == 20:
                print forces[natoms*3:]
                print np.amax(np.absolute(np.subtract(self.qcell[0:natoms*3], self.gmc.xold[0:natoms*3])))
            print 'gradient', g[natoms*3:]
         
        self.qcell[0:natoms*3] = self.beads.q
        self.qcell[natoms*3:] = self.gmc.strain*self.jacobian
        #self.qcell = self.gmc.strain*self.jacobian
        alam = 0.5
        self.qcell = self.qcell + alam*g
        self.beads.q = self.qcell[0:natoms*3]
        eps_vec = self.qcell[natoms*3:]/self.jacobian
        #eps_vec = self.qcell/self.jacobian
        #unit = np.eye(3, dtype=float)
        #eps = np.reshape(eps_vec, (3,3)) + unit
        eps = np.array([[1.0 + eps_vec[0], 0.5 * eps_vec[5], 0.5 * eps_vec[4]],
                        [0.5 * eps_vec[5], 1.0 + eps_vec[1], 0.5 * eps_vec[3]],
                        [0.5 * eps_vec[4], 0.5 * eps_vec[3], 1.0 + eps_vec[2]]])
        self.cell.h = np.dot(self.gmc.oldcell,eps)
        self.gmc.strain = eps_vec
        
        self.gmc.xold[:] = self.qcell
        self.gmc.oldcell = self.cell.h.copy()
        
        print 'strain', self.gmc.strain
        print 'newcell', self.cell.h
        
        self.qtime += time.time()
        print self.forces.pot
        
        if step ==20:
            softexit.trigger('Teststep')


class CellOptimizer(DummyOptimizer):

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
                
        if step == 0:   # or np.sqrt(np.dot(self.gm.d, self.gm.d)) == 0.0: this part for restarting at claimed minimum (optional)
            info(" @GEOP: Initializing BFGS", verbosity.debug)
            
            natoms = self.beads.natoms
            # new force vector with forces and transformed stress tensor
            stress = self.forces.vir.copy()/self.cell.V*self.norm_stress
            # gradient direction
            self.gmc.d = stress.flatten()/ np.sqrt(np.dot(stress.flatten(), stress.flatten()))
            # store actual position to previous position, initial strain is zero
            self.gmc.xold = np.zeros(9, float)
            self.gmc.oldcell = self.cell.h.copy()
            self.gmc.strain = np.zeros(9, float)
            
        natoms = self.beads.natoms
        # Current energy and forces
        u0 = self.forces.pot.copy()
        # new force vector with forces and transformed stress tensor
        stress = self.forces.vir.copy()/self.cell.V*self.norm_stress
        #forces = self.forces.vir.flatten()*self.norm_stress/self.cell.V
        du0 = - stress.flatten()
        print du0
        print 'Strain', self.gmc.strain

        # Store previous forces
        self.old_f[:] = stress.flatten()
        #print 'before BFGS', self.gmc.xold[0], self.gmc.oldcell[0,0]
        
        self.qcell = self.gmc.strain*self.jacobian
        
        # Do one iteration of BFGS, return new point, function value,
        # move direction, and current Hessian to use for next iteration
        self.qcell, fx, self.gmc.d, self.invhessian = BFGS(self.qcell,
                self.gmc.d, self.gmc, fdf0=(u0, du0), invhessian=self.invhessian,
                big_step=self.big_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"])
        
        eps_vec = self.qcell/self.jacobian
        eps = np.reshape(eps_vec, (3,3))
        #eps = np.array([[1.0 + eps_vec[0], 0.5 * eps_vec[1], 0.5 * eps_vec[2]],
        #                [0.5 * eps_vec[1], 1.0 + eps_vec[4], 0.5 * eps_vec[5]],
         #               [0.5 * eps_vec[2], 0.5 * eps_vec[5], 1.0 + eps_vec[8]]])
        unit = np.eye(3, dtype=float)
        self.cell.h = np.dot(self.gmc.oldcell, unit + eps)
        #self.cell.h = np.dot(self.gmc.oldcell, eps)
        print 'Cell', self.cell.h
        self.gmc.strain = eps_vec
        
        print 'strain', self.gmc.strain
        
        x = np.amax(np.absolute(np.subtract(self.cell.h,self.gmc.oldcell)))
        print 'x', x
        # Store old position
        self.gmc.xold[:] = self.qcell
        self.gmc.oldcell = self.cell.h.copy()
                
        stress = self.forces.vir.copy()/self.cell.V*self.norm_stress
        
        # Exit simulation step
        self.exitstep2(fx, u0, x, stress)
'''        
class BFGSCellOptimizer(DummyOptimizer):
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
                
        if step == 0:   # or np.sqrt(np.dot(self.gm.d, self.gm.d)) == 0.0: this part for restarting at claimed minimum (optional)
            info(" @GEOP: Initializing BFGS", verbosity.debug)
            
            #number of atoms
            natoms = self.beads.natoms
            # new force vector with forces and transformed stress tensor
            forces = np.zeros((natoms+2)*3, float)
            forces[0:natoms*3] = self.forces.f
            stress = self.forces.vir.flatten()/self.cell.V
            stressmat = np.array([stress[0], stress[4], stress[8],
                        stress[5], stress[2], stress[1]])
            forces[natoms*3:] = stressmat*self.norm_stress
            #print 'forces', forces[0:natoms*3]*51.422065
            #print 'stress', stressmat* 183.63153
            #print 'pos', self.beads.q*0.52917721
            #print 'cell', self.cell.h*0.52917721
            #softexit.trigger('aus')
            # gradient direction
            self.gmc.d = forces/ np.sqrt(np.dot(forces.flatten(), forces.flatten()))
            # store actual position to previous position, initial strain is zero
            self.gmc.xold = np.zeros((natoms+2)*3, float)
            self.gmc.xold[0:natoms*3] = self.beads.q.copy()
            self.gmc.oldcell = self.cell.h.copy()
            self.gmc.strain = np.zeros(6, float)
        
        #number of atoms   
        natoms = self.beads.natoms
        # Current energy and forces
        u0 = self.forces.pot.copy()
        # new force vector with forces and transformed stress tensor
        forces = np.zeros((natoms+2)*3,float)
        forces[0:natoms*3] = self.forces.f
        stress = self.forces.vir.flatten()/self.cell.V
        stressmat = np.array([stress[0], stress[4], stress[8],
                        stress[5], stress[2], stress[1]])
        forces[natoms*3:] = stressmat*self.norm_stress
        #forces[natoms*3:] = self.forces.vir.flatten()*self.norm_stress/self.cell.V
        du0 = - forces

        # Store previous forces
        self.old_f[:] = forces
        
        #new vector containing atomic positions and transformed strain
        self.qcell[0:natoms*3] = self.beads.q
        self.qcell[natoms*3:] = 0.0#self.gmc.strain*self.jacobian
        
        # Do one iteration of BFGS, return new point, function value,
        # move direction, and current Hessian to use for next iteration
        self.qcell, fx, self.gmc.d, self.invhessian, forces = BFGS(self.qcell,
                self.gmc.d, self.gmc, fdf0=(u0, du0), invhessian=self.invhessian,
                big_step=self.big_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"])
        
        # x = current position - previous position; use for exit tolerance
        x = np.amax(np.absolute(np.subtract(self.qcell, self.gmc.xold)))
        
        # updating positions and cell by transforming back
        self.beads.q[-1,:] = self.qcell[0:natoms*3]
        eps_vec = self.qcell[natoms*3:]/self.jacobian
        #eps = np.reshape(eps_vec, (3,3))
        eps = np.array([[1.0 + eps_vec[0], 0.5 * eps_vec[5], 0.5 * eps_vec[4]],
                        [0.5 * eps_vec[5], 1.0 + eps_vec[1], 0.5 * eps_vec[3]],
                        [0.5 * eps_vec[4], 0.5 * eps_vec[3], 1.0 + eps_vec[2]]])
        #unit = np.eye(3, dtype=float)
        #self.cell.h = np.dot(self.gmc.oldcell, unit + eps)
        self.cell.h = np.dot(self.gmc.oldcell, eps)
        self.gmc.strain = eps_vec
        print self.gmc.strain
        
        #M = np.linalg.solve(self.gmc.oldcell, self.cell.h)
        #print 'M', M
        #qresh=np.reshape(self.qcell[0:natoms*3],(natoms,3))
        #for i in range(natoms):
         #   qresh[i,:]=np.dot(qresh[i,:],M)
        #self.beads.q[-1,:] = qresh.flatten()

        # Store old position
        self.gmc.xold = self.qcell.copy()
        self.gmc.oldcell = self.cell.h.copy()
        
        #forces = np.zeros((natoms+2)*3,float)
        #forces[0:natoms*3] = self.forces.f
        #stress = self.forces.vir.flatten()/self.cell.V
        #stressmat = np.array([[stress[0], stress[4], stress[8]],
         #               [stress[5], stress[2], stress[1]]])
        #forces[natoms*3:] = stressmat.flatten()*self.norm_stress
        
        if step == 20:
            softexit.trigger('step2')
        # Exit simulation step
        self.exitstep2(fx, u0, x, forces)

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
         
        if step == 0:   # or np.sqrt(np.dot(self.gm.d, self.gm.d)) == 0.0: this part for restarting at claimed minimum (optional)
            info(" @GEOP: Initializing BFGS", verbosity.debug)
            self.gm.d = depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
            # store actual position to previous position
            self.gm.xold = self.beads.q.copy()
        
        # Current energy and forces
        u0 = self.forces.pot.copy()
        du0 = - self.forces.f

        # Store previous forces
        self.old_f[:] = self.forces.f

        # Do one iteration of BFGS, return new point, function value,
        # move direction, and current Hessian to use for next iteration
        self.beads.q, fx, self.gm.d, self.invhessian = BFGS(self.beads.q,
                self.gm.d, self.gm, fdf0=(u0, du0), invhessian=self.invhessian,
                big_step=self.big_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"])
                
        # x = current position - previous position; use for exit tolerance
        x = np.amax(np.absolute(np.subtract(self.beads.q, self.gm.xold)))
        
        
        # Store old position
        self.gm.xold[:] = self.beads.q
        
        # Exit simulation step
        self.exitstep(fx, u0, x)

class LBFGSOptimizer(DummyOptimizer):
    """ L-BFGS Minimization """
    
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
        
        # Initialize approximate Hessian inverse to the identity and direction
        # to the steepest descent direction
        if step == 0:   # or np.sqrt(np.dot(self.gm.d, self.gm.d)) == 0.0: <-- this part for restarting at claimed minimum (optional)
            info(" @GEOP: Initializing L-BFGS", verbosity.debug)
            self.gm.d = depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
            # store actual position to previous position
            self.gm.xold = self.beads.q.copy()
            # Initialize lists of previous positions and gradient
            self.qlist = np.zeros((self.corrections, len(self.beads.q.flatten())))
            self.glist = np.zeros((self.corrections, len(self.beads.q.flatten())))

        # Current energy and force
        u0, du0 = (self.forces.pot.copy(), - self.forces.f)

        # Store previous forces
        self.old_f[:] = self.forces.f

        # Do one iteration of L-BFGS, return new point, function value,
        # move direction, and current Hessian to use for next iteration
        self.beads.q, fx, self.gm.d, self.qlist, self.glist = L_BFGS(self.beads.q,
                self.gm.d, self.gm, self.qlist, self.glist,
                fdf0=(u0, du0), big_step=self.big_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"],
                m=self.corrections, k=step)

        info(" @GEOP: Updated position list", verbosity.debug)
        info(" @GEOP: Updated gradient list", verbosity.debug)

        # x = current position - old position. Used for convergence tolerance
        x = np.amax(np.absolute(np.subtract(self.beads.q, self.gm.xold)))
        
        # Store old position
        self.gm.xold[:] = self.beads.q
        
        # Exit simulation step
        self.exitstep(fx, u0, x)
         
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

        gradf1 = dq1 = depstrip(self.forces.f)

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
        
        self.lm.set_dir(depstrip(self.beads.q), dq1_unit)
        
        # Reuse initial value since we have energy and forces already
        u0, du0 = (self.forces.pot.copy(), np.dot(depstrip(self.forces.f.flatten()), dq1_unit.flatten()))

        # Do one SD iteration; return positions and energy
        (x, fx) = min_brent(self.lm, fdf0=(u0, du0), x0=0.0,
                    tol=self.ls_options["tolerance"],
                    itmax=self.ls_options["iter"], init_step=self.ls_options["step"])
        # Automatically adapt the search step for the next iteration.
        # Relaxes better with very small step --> multiply by factor of 0.1 or 0.01
        self.ls_options["step"] = 0.1 * x * self.ls_options["adaptive"] + (1 - self.ls_options["adaptive"]) * self.ls_options["step"]

        self.beads.q += dq1_unit * x
        
        # Exit simulation step
        self.exitstep(fx, u0, x)

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

        if step == 0:
            gradf1 = dq1 = depstrip(self.forces.f)

            # Move direction for 1st conjugate gradient step
            dq1_unit = dq1 / np.sqrt(np.dot(gradf1.flatten(), gradf1.flatten()))
            info(" @GEOP: Determined SD direction", verbosity.debug)
    
        else:
        
            gradf0 = self.old_f
            dq0 = self.old_d
            gradf1 = depstrip(self.forces.f)
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

        self.lm.set_dir(depstrip(self.beads.q), dq1_unit)

        # Reuse initial value since we have energy and forces already
        u0, du0 = (self.forces.pot.copy(), np.dot(depstrip(self.forces.f.flatten()), dq1_unit.flatten()))

        # Do one CG iteration; return positions and energy
        (x, fx) = min_brent(self.lm, fdf0=(u0, du0), x0=0.0,
                    tol=self.ls_options["tolerance"],
                    itmax=self.ls_options["iter"], init_step=self.ls_options["step"])

        # Automatically adapt the search step for the next iteration.
        # Relaxes better with very small step --> multiply by factor of 0.1 or 0.01
        self.ls_options["step"] = 0.1 * x * self.ls_options["adaptive"] + (1 -     self.ls_options["adaptive"]) * self.ls_options["step"]

        self.beads.q += dq1_unit * x
        
        # Exit simulation step
        self.exitstep(fx, u0, x)
