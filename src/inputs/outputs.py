"""Deals with creating the output objects.

Classes:
   InputOutputs: Creates a list of all the output objects.
   InputProperties: Deals with property output.
   InputTrajectory: Deals with trajectory output.
   InputCheckpoint: Deals with restart file output.
"""

from utils.depend import *
from utils.inputvalue import *
from copy import copy
import engine.outputs
import numpy as np
from engine.properties import getkey

__all__=['InputOutputs', 'InputProperties', 'InputTrajectory',
         'InputCheckpoint']

class InputProperties(InputArray):
   """Simple input class to describe output for properties.

   Storage class for PropertyOutput.

   Attributes:
      filename: The name of the file to output to.
      stride: The number of steps that should be taken between outputting the
         data to file.
   """

   default_help = """This class deals with the output of one property. """
   default_label = "PROPERTIES"

   attribs=copy(InputArray.attribs)
   attribs["filename"]=(InputAttribute,{ "dtype" : str, "default": "out"} )
   attribs["stride"]=(InputAttribute,{ "dtype" : int, "default": 1 } )

   def __init__(self, help=None,  default=None, dtype=None, dimension=None):
      """Initializes InputProperties.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputProperties,self).__init__(help=help, default=default, dtype=str, dimension=dimension)

   def fetch(self):
      """Returns a PropertyOutput object."""

      return engine.outputs.PropertyOutput(self.filename.fetch(), self.stride.fetch(), super(InputProperties,self).fetch())

   def store(self, prop):
      """Stores a PropertyOutput object."""

      super(InputProperties,self).store(prop.outlist)
      self.stride.store(prop.stride)
      self.filename.store(prop.filename)


class InputTrajectory(InputValue):
   """Simple input class to describe output for trajectories.

   Storage class for TrajectoryOutput.

   Attributes:
      filename: The (base) name of the file to output to.
      stride: The number of steps that should be taken between outputting the
         data to file.
      format: The format of the trajectory output file.
   """

   default_help = """This class defines how one trajectory file should be output. """
   default_label = "TRAJECTORY"

   attribs=copy(InputValue.attribs)
   attribs["filename"]=(InputAttribute,{ "dtype" : str, "default": "traj"} )
   attribs["stride"]=(InputAttribute,{ "dtype" : int, "default": 1 } )
   attribs["format"]=(InputAttribute,{ "dtype" : str, "default": "xyz" } )

   def __init__(self, help=None,  default=None, dtype=None, dimension=None):
      """Initializes InputTrajectory.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputTrajectory,self).__init__(help=help, default=default, dtype=str, dimension=dimension)

   def fetch(self):
      """Returns a TrajectoryOutput object."""

      return engine.outputs.TrajectoryOutput(self.filename.fetch(), self.stride.fetch(), super(InputTrajectory,self).fetch(),self.format.fetch())

   def store(self, traj):
      """Stores a PropertyOutput object."""

      super(InputTrajectory,self).store(traj.what)
      self.stride.store(traj.stride)
      self.filename.store(traj.filename)
      self.format.store(traj.format)


class InputCheckpoint(InputValue):
   """Simple input class to describe output for properties.

   Storage class for CheckpointOutput.

   Attributes:
      filename: The (base) name of the file to output to.
      stride: The number of steps that should be taken between outputting the
         data to file.
      overwrite: whether checkpoints should be overwritten, or multiple
         files output.
   """

   default_help = """This class defines how a checkpoint file should be output. """
   default_label = "CHECKPOINT"

   attribs=copy(InputValue.attribs)
   attribs["filename"]=(InputAttribute,{ "dtype" : str, "default": "restart"} )
   attribs["stride"]=(InputAttribute,{ "dtype" : int, "default": 1 } )
   attribs["overwrite"]=(InputAttribute,{ "dtype" : bool, "default": True } )

   def __init__(self, help=None,  default=None, dtype=None, dimension=None):
      """Initializes InputCheckpoint.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputCheckpoint,self).__init__(help=help, default=default, dtype=int, dimension=dimension)

   def fetch(self):
      """Returns a CheckpointOutput object."""

      step=super(InputCheckpoint,self).fetch()
      return engine.outputs.CheckpointOutput(self.filename.fetch(), self.stride.fetch(), self.overwrite.fetch(), step=step )

   def parse(self, xml=None, text=""):
      """Overwrites the standard parse function so that we can specify this tag
      in the input without any data.

      We can use the syntax <checkpoint /> to do this

      Args:
         xml: An xml node containing all the data for the parent tag.
         text: The data to read the data from. Will be None if we have not
            specified any data.
      """

      # just a quick hack to allow an empty element
      try:
         super(InputCheckpoint,self).parse(xml,text)
      except: #TODO make this except a specific exception, not every one
         self.value=0  #This could hide actual errors, at least in theory.

   def store(self, chk):
      """Stores a CheckpointOutput object."""

      super(InputCheckpoint,self).store(chk.step)
      self.stride.store(chk.stride)
      self.filename.store(chk.filename)
      self.overwrite.store(chk.overwrite)


class InputOutputs(Input):
   """ List of outputs input class.

   An example of a dynamic input class: a variable number of tags might be
   present, corresponding to different output requests. This allows for
   instance to print multiple property outputs, with different content
   and/or output frequency.

   Attributes:
      prefix: A string that will be appended to all output files from this
         simulation.
      extra: A list of all the output objects.
   """

   attribs = { "prefix" : ( InputAttribute, { "dtype" : str,
                                          "default"  : "wrap-pi",
                                          "help"     : "A string that will be the pre-pended to each output file name." })
             }

   dynamic = {  "properties" : (InputProperties, { "help" : "Each of the <properties> tags specify how to create a file in which one or more properties are written, one line per frame. " } ),
               "trajectory" : (InputTrajectory, { "help" : "Each of the <trajectory> tags specify how to create a trajectory file, containing a list of per-atom-coordinate properties. " } ),
               "checkpoint" : (InputCheckpoint, { "help" : "Each of the <checkpoint> tags specify how to create a checkpoint file, which can be used to restart a simulation. " } ),
            }

   default_help = """This class defines how properties, trajectories and checkpoints should be output during the simulation.
    May contain zero, one or many instances of <properties>, <trajectory> or <checkpoint> tags, each giving instructions on how
    one output file should be created and managed. """
   default_label = "OUTPUTS"

   @classmethod
   def make_default(cls):
      """Used to make the default value of the outputs class for use when no
      output is specified.

      Needed since this is a fairly complicated default, with many mutable
      objects, and the default has to be generated by a function that does not
      use any mutable objects as arguments.
      """

      return [ engine.outputs.PropertyOutput("wrap-pi.md", 10, [ "time", "step", "conserved", "temperature", "potential", "kinetic_cv" ] ),
               engine.outputs.TrajectoryOutput("wrap-pi.pos", 100, "positions", "xyz"),
               engine.outputs.CheckpointOutput("wrap-pi.checkpoint",1000,overwrite=True)]

   def fetch(self):
      """Returns a list of the output objects included in this dynamic
      container.

      Returns:
         A list of tuples, with each tuple being of the form ('type', 'object')
         where 'type' is the type of output object and 'object' is a particular
         object of that type.
      """

      super(InputOutputs, self).fetch()
      outlist = [ p.fetch() for (n, p) in self.extra ]
      prefix = self.prefix.fetch()
      if not prefix == "":
         for p in outlist:
            p.filename = prefix + "." + p.filename

      return outlist

   def store(self, plist):
      """ Stores a list of the output objects, creating a sequence of
      dynamic containers.

      Args:
         plist: A list of tuples, with each tuple being of the form
            ('type', 'object') where 'type' is the type of forcefield and
            'object' is a particular object of that type.
      """

      super(InputOutputs, self).store()
      self.extra = []

      self.prefix.store("")
      for el in plist:
         if (isinstance(el, engine.outputs.PropertyOutput)):
            ip = InputProperties()
            ip.store(el)
            self.extra.append(("properties", ip))
         elif (isinstance(el, engine.outputs.TrajectoryOutput)):
            ip = InputTrajectory()
            ip.store(el)
            self.extra.append(("trajectory", ip))
         elif (isinstance(el, engine.outputs.CheckpointOutput)):
            ip = InputCheckpoint()
            ip.store(el)
            self.extra.append(("checkpoint", ip))
