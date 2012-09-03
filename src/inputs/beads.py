"""Deals with creating the beads class.

Classes:
   InputStartBeads: Deals with starting a simulation from a different
      number of beads.
   InputBeads: Deals with creating the Beads object from a file, and
      writing the checkpoints.
"""

import numpy as np
import math
from engine.beads import *
from engine.atoms import Atoms
from utils.inputvalue import *
import utils.io.io_pdb
from utils.depend import *
from utils.units import *
from inputs.atoms import *

__all__ = ['InputBeads']

class InputBeads(Input):
   """Beads input class.

   Handles generating the appropriate beads class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      nbeads: An optional integer giving the number of beads. Defaults to 0.
      natoms: An optional integer giving the number of atoms. Defaults to 0.


   Fields:
      start_centroid: An atoms object to initialize the centroid postions from.
      start_beads: A beads object to initialize the normal mode
         coordinates from.
      q: An optional array giving the bead positions. Defaults to an empty
         array with no elements.
      p: An optional array giving the bead momenta. Defaults to an empty
         array with no elements.
      m: An optional array giving the bead masses. Defaults to an empty array
         with no elements.
      names: An optional array giving the bead names. Defaults to an empty
         array with no elements.
   """

   attribs = {  "natoms"    : (InputAttribute, {"dtype"     : int,  "default"   : 0,
                                            "help"      : "The number of atoms."}),
                "nbeads"    : (InputAttribute, {"dtype"     : int,  "default"   : 0,
                                            "help"      : "The number of beads."})
            }
   fields={ "q"         : (InputArray, {"dtype"     : float,
                                        "default"   : input_default(factory=np.zeros, args = (0,)),
                                        "help"      : "The positions of the beads. In an array of size [nbeads, 3*natoms].",
                                        "dimension" : "length"}),
            "p"         : (InputArray, {"dtype"     : float,
                                        "default"   : input_default(factory=np.zeros, args = (0,)),
                                        "help"      : "The momenta of the beads. In an array of size [nbeads, 3*natoms].",
                                        "dimension" : "momentum"}),
            "m"         : (InputArray, {"dtype"     : float,
                                        "default"   : input_default(factory=np.zeros, args = (0,)),
                                        "help"      : "The masses of the atoms, in the format [m1, m2, ... ].",
                                        "dimension" : "mass"}),
            "names"     : (InputArray, {"dtype"     : str,
                                        "default"   : input_default(factory=np.zeros, args=(0,), kwargs={'dtype': np.dtype('|S6')}),
                                        "help"      : "The names of the atoms, in the format [name1, name2, ... ]."})  }

   default_help = "Describes the configurations of atoms in a path integral simulations."
   default_label = "BEADS"


   def store(self, beads):
      """Takes a Beads instance and stores a minimal representation of it.

      Args:
         beads: A Beads object from which to initialise from.
      """

      super(InputBeads,self).store()
      self.natoms.store(beads.natoms)
      self.nbeads.store(beads.nbeads)

      self.q.store(depstrip(beads.q))
      self.p.store(depstrip(beads.p))
      self.m.store(depstrip(beads.m))
      self.names.store(depstrip(beads.names))

   def fetch(self):
      """Creates a beads object.

      Returns:
         A beads object of the appropriate type and with the appropriate
         properties given the attributes of the InputBeads object.
      """

      super(InputBeads,self).fetch()
      beads = Beads(self.natoms.fetch(),self.nbeads.fetch())

      # tries to fill up with as much data as available and valid
      # TODO should print warnings (or just abort?) if the arrays are not empty and size mismatch
      q=self.q.fetch();      
      if (q.shape == (beads.nbeads,3*beads.natoms)) : beads.q=q
      elif len(q) != 0 : raise ValueError("Array shape mismatches for q in <beads> input.")
      p=self.p.fetch();
      if (p.shape == (beads.nbeads,3*beads.natoms)) : beads.p=p
      elif len(p) != 0 : raise ValueError("Array shape mismatches for p in <beads> input.")
      m=self.m.fetch();
      if (m.shape == (beads.natoms,)) : beads.m=m
      elif len(m) != 0 : raise ValueError("Array shape mismatches for m in <beads> input.")
      n=self.names.fetch();
      if (n.shape == (beads.natoms,)) : beads.names=n
      elif len(n) != 0 : raise ValueError("Array shape mismatches for names in <beads> input.")

      return beads
