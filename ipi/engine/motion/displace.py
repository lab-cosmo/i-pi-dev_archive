# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import time
import numpy as np

from ipi.engine.motion import Motion
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import unit_to_internal


__all__ = ['Displace']


class Displace(Motion):
    """Calculator object that just displaces all non-fixed atoms by a specific amount for each step

    Attributes:
        dv: the rigid displacement vector

    Depend objects:
        None really meaningful.
    """

    def __init__(self, fixcom=False, fixatoms=None, displacement=None):
        """Initialises Replay.

        Args:
           temp: The system temperature.
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
           displacement: displacement vector of coordinates
        """
        super(Displace, self).__init__(fixcom=fixcom, fixatoms=fixatoms)
        if displacement is None:
            raise ValueError("Must provide a displacement vector")
        self.displacement = displacement
        if self.displacement.size != 3:
            raise ValueError("Displacement vector must be 3 values")

    def bind(self, ens, beads, nm, cell, bforce, prng):
        super(Displace, self).bind(ens, beads, nm, cell, bforce, prng)

    def step(self, step=None):
        """Step the displacement."""
        if step is None:
            softexit.trigger("Unknown step count. Exiting simulation")
            return

        for bead in self.beads:
            # To get correct displacement of all atoms
            n = bead.q.size / 3
            bead.q += np.tile(self.displacement, n)
