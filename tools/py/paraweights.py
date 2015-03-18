#!/usr/bin/python
""" parasort.py

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

Post-processes the output of a parallel-tempering simulation and
re-orders the outputs so that they correspond to the different
temperatures ensembles rather than to the time series of one of
the replicas exchanging temperature over time. 
Given a target temperature for re-weighing, it will also print out 
relative weights for the different trajectories based on 
Ceriotti, Brain, Riordan, Manolopoulos, Proc. Royal Soc. A 468 (2011)

It should be run in the same dyrectory as where i-pi was (or is being)
run, and simply fetches all information from the simulation input file.
Will create a series of PTindex_* files, each corresponding to the
data for replica 'index'.

Syntax:
   parasort.py inputfile.xml [prefix] [temperature(K)]
"""

import sys, re
import numpy as np
from ipi.engine.simulation import Simulation
from ipi.engine.outputs import *
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.io_xml import *
from ipi.utils.units import unit_to_internal
from ipi.utils.mathtools import logsumlog


def main(inputfile, prefix="PTW-", ttemp="300.0", skip="5000"):
   txtemp = ttemp
   ttemp = unit_to_internal("energy","kelvin",float(ttemp))
   skip = int(skip)
   
   # opens & parses the input file
   ifile = open(inputfile,"r")
   xmlrestart = xml_parse_file(ifile) # Parses the file.
   ifile.close()

   isimul = InputSimulation()
   isimul.parse(xmlrestart.fields[0][1])

   simul = isimul.fetch()

   if simul.mode != "paratemp":
      raise ValueError("Simulation does not look like a parallel tempering one.")

   # reconstructs the list of the property and trajectory files that have been output
   # and that should be re-ordered
   lprop = [ ] # list of property files
   ltraj = [ ] # list of trajectory files
   tlist = simul.paratemp.temp_list
   nsys = len(simul.syslist)
   for o in simul.outtemplate:
      if type(o) is CheckpointOutput:   # properties and trajectories are output per system
         pass
      elif type(o) is PropertyOutput:
         nprop =  []
         isys=0
         for s in simul.syslist:   # create multiple copies
            if s.prefix != "":
               filename = s.prefix+"_"+o.filename
            else: filename=o.filename
            ofilename = prefix+ ( ("%0" + str(int(1 + np.floor(np.log(len(simul.syslist))/np.log(10)))) + "d") % (isys) )+"_"+o.filename
            nprop.append( { "filename" : filename, "ofilename" : ofilename, "stride": o.stride,
                           "ifile" : open(filename, "r"), "ofile" : None
             } )
            isys+=1
         lprop.append(nprop)         

   ptfile=open("PARATEMP", "r")

   # these are variables used to compute the weighting factors 
   tprops=[]
   vfields=[]
   vunits=[]
   nw = []
   tw = [] 
   tw2 = []
   
   repot = re.compile(' ([0-9]*) *--> potential')
   reunit = re.compile('{(.*)}')
   
 
   # now reads files one frame at a time, and re-direct output to the appropriate location
   irep = np.zeros(nsys,int)
   while True:
      # reads one line from PARATEMP index file
      line=ptfile.readline()
      line = line.split()
      if len(line) == 0: break

      step = int(line[0])
      irep[:] = line[1:]

      try:
         wk = 0
         for prop in lprop:
            for isys in range(nsys):
               sprop = prop[isys]
               if step % sprop["stride"] == 0: # property transfer
                  iline = sprop["ifile"].readline()
                  if len(iline)==0: raise EOFError
                  while iline[0] == "#":  # fast forward if line is a comment 
                     # checks if we have one single file with potential energies
                     rm=repot.search(iline)
                     if not (rm is None) and not (prop in tprops):
                        tprops.append(prop)
                        for p in prop:
                           p["ofile"] = open(p["ofilename"],"w")
                           p["ofile"].write("# column   1     --> ptlogweight: ln of re-weighing factor with target temperature %s K\n" % (txtemp) )
                        vfields.append(int(rm.group(1)))
                        nw.append(np.zeros(nsys))
                        tw.append(np.zeros(nsys))
                        tw2.append(np.zeros(nsys))
                        rm=reunit.search(iline)
                        if rm: vunits.append(rm.group(1))
                        else: vunits.append("atomic_unit")                              
                     iline = sprop["ifile"].readline()                     
                  if prop in tprops: # do temperature weighing    
                     pot = unit_to_internal("energy", vunits[wk],float(iline.split()[vfields[wk]]))
                     ir = irep[isys]
                     temp = tlist[ir]
                     lw = pot*(1/temp-1/ttemp)
                     if step > skip: # computes trajectory weights avoiding the initial - possibly insane - values
                        if nw[wk][ir] ==0 :
                           tw[wk][ir] = lw
                           tw2[wk][ir] = lw
                        else:
                           tw[wk][ir] = logsumlog((tw[wk][ir],1),(lw,1))[0]
                           tw2[wk][ir] = logsumlog((tw2[wk][ir],1),(2*lw,1))[0]                     
                        nw[wk][ir] += 1
                     prop[ir]["ofile"].write("%15.7e\n" %(lw))
                     if isys == nsys-1: wk+=1
      except EOFError:
         # print out trajectory weights based on PRSA 2011, assuming that observables and weights are weakly correlated
         wk=0
         fpw = open(prefix+"WEIGHTS","w")
         for prop in lprop:
            if prop in tprops:
               for ir in range(nsys):                  
                  fpw.write("%s   %15.7e \n" % (prop[ir]["ofilename"], 1.0/(np.exp(tw2[wk][ir]-2*tw[wk][ir])*nw[wk][ir])  )  )
               wk+=1
         break

if __name__ == '__main__':
   main(*sys.argv[1:])
