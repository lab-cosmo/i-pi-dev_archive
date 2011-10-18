import numpy, math, sys
import cell 

def print_pdb(atoms, ncell, filedesc = sys.stdout):
   """Takes the system and gives pdb formatted output for the unit cell and the
      atomic positions """

   a, b, c, alpha, beta, gamma = cell.h2abc(ncell.h)
   alpha *= 180.0/math.pi #radian to degree conversion
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi
   
   z = 1 #number of polymeric chains in a unit cell. I can't decide if 1 or 0 is more sensible for this...

   filedesc.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a, b, c, alpha, beta, gamma, " P 1        ", z))

   for i in range(0,len(atoms)): 
      filedesc.write("ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n" % (i+1, atoms[i].name,' ','  1',' ',1,' ',atoms[i].q[0],atoms[i].q[1],atoms[i].q[2],0.0,0.0,'  ',0))

   filedesc.write("END\n")





def read_pdb(filedesc):
   """Takes a pdb-style file and creates a system with the appropriate unit
      cell and atom positions"""

   header = filedesc.readline()
   a = float(header[6:15]);      b = float(header[15:24]);    c = float(header[24:33]);
   alpha = float(header[33:40]); beta = float(header[40:47]); gamma = float(header[47:54]);
   alpha *= math.pi/180.0;       beta *= math.pi/180.0;       gamma *= math.pi/180.0
   cell = numpy.array([a, b, c, alpha, beta, gamma])

   atoms = []
   natoms = 0
   body = filedesc.readline()
   while body != '':
      natoms += 1
      name = body[12:16]
      x = float(body[31:39])
      y = float(body[39:47])
      z = float(body[47:55])
      pos = numpy.array([x, y, z])
      atoms.append([name, pos]) 
      body = filedesc.readline()
   return atoms, cell, natoms

def xml(system, namedpipe):
   tab = "   "
   namedpipe.write("<?xml version=\"1.0\"?>\n")
   namedpipe.write("<system>\n")
   namedpipe.write(tab + "<natoms>" + str(system.natoms) + "</natoms>\n")
   namedpipe.write(tab + "<temp>" + str(system.temp) + "</temp>\n")
   namedpipe.write(tab + "<cutoff>" + str(system.cell.cutoff) + "</cutoff>\n")
   namedpipe.write(tab + "<atoms>\n")

   for i in range(system.natoms):
      atom_q = system.atoms[i].q
      namedpipe.write(tab + tab + "<atom>\n")
      namedpipe.write(tab + tab + tab + "<index>" + str(i) + "</index>\n")
      namedpipe.write(tab + tab + tab + "<x>" + str(atom_q[0]) + "</x>\n")
      namedpipe.write(tab + tab + tab + "<y>" + str(atom_q[1]) + "</y>\n")
      namedpipe.write(tab + tab + tab + "<z>" + str(atom_q[2]) + "</z>\n")
      namedpipe.write(tab + tab + tab + "<mass>" + str(system.atoms[i].mass) + "</mass>\n")
      namedpipe.write(tab + tab + "</atom>\n")
   namedpipe.write(tab + "</atoms>\n")

   h = system.cell.h
   namedpipe.write(tab + "<cell>\n")

   namedpipe.write(tab + tab + "<h1>\n")
   namedpipe.write(tab + tab + tab + "<x>" + str(h[0,0]) + "</x>\n")
   namedpipe.write(tab + tab + tab + "<y>" + str(0.0) + "</y>\n")
   namedpipe.write(tab + tab + tab + "<z>" + str(0.0) + "</z>\n")
   namedpipe.write(tab + tab + "</h1>\n")

   namedpipe.write(tab + tab + "<h2>\n")
   namedpipe.write(tab + tab + tab + "<x>" + str(h[0,1]) + "</x>\n")
   namedpipe.write(tab + tab + tab + "<y>" + str(h[1,1]) + "</y>\n")
   namedpipe.write(tab + tab + tab + "<z>" + str(0.0) + "</z>\n")
   namedpipe.write(tab + tab + "</h2>\n")
   
   namedpipe.write(tab + tab + "<h3>\n")
   namedpipe.write(tab + tab + tab + "<x>" + str(h[0,2]) + "</x>\n")
   namedpipe.write(tab + tab + tab + "<y>" + str(h[1,2]) + "</y>\n")
   namedpipe.write(tab + tab + tab + "<z>" + str(h[2,2]) + "</z>\n")
   namedpipe.write(tab + tab + "</h3>\n")

   namedpipe.write(tab + "</cell>\n")

   namedpipe.write("</system>\n") 
