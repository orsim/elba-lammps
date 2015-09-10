#!/usr/bin/env python

# Script: trj2pdb-extractWater.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Reads a LAMMPS ".trj" trajectory (dump) file and converts it
#          into a VMD-compatible file "trajectory.pdb". Dedicated 'state'
#          files can be found in 'elba-lammps/viz/'.
# Syntax: trj2pdb-extractWater.py inputFile
# Example: trj2pdb-extractWater.py dump.trj
# Notes: Point-dipoles in LAMMPS ".trj" files are described by mass center
#        coordinates and orientations (x-, y-, and z- projections of the
#        dipole vectors). Since VMD does not allow visualization of point-
#        dipoles, this script uses the orientation information to convert
#        every dipolar atom into two atoms representing the "+" and "-"
#        ends of the original dipole. Such two atoms are spaced 1 Angstrom
#        apart. 
# Preprocessing: sometimes it may be useful to concatenate individual 
#                LAMMPS trajectory "snapshots": cat dump.*.trj > dump.trj 
# Reference: http://www.ks.uiuc.edu/Research/vmd

import sys, string, linecache
from math import sqrt

if len(sys.argv) != 2:
  print "Syntax: trj2pdb-extractWater.py inputFile"
  sys.exit()

inFileName = sys.argv[1]
inFile = open(inFileName, "r")

print "Processing file %s ..." % inFileName
line = linecache.getline(inFileName, 4)
words = string.split(line)
nAtoms = int(words[0])
print "Number of \'atoms\': %d" % nAtoms

outFile = open("elbawater.pdb", "w")
lines = inFile.readlines()

SCALE=0.05; # [Angstrom] scaling factor
nLine=0; # line counter
nCount=0; # counter
nExtra=0; # extra atoms used to represent dipoles
watType = 9999999 # initilize water type

print "Scanning coordinates and orientations..."
for line in lines:
    nLine=nLine+1
    words = string.split(line)
    if len(words) == 2:
        if words[1]=="TIMESTEP":
            if nLine!=1:
                outFile.write('ENDMDL\n')
                nExtra=0;
            outFile.write('MODEL\n')
            outFile.write('CRYST1')
        else:
            nCount=nCount+1
            negCoord=float(words[0])
            plusCoord=float(words[1])
            boxLength=plusCoord-negCoord
            outFile.write('%9.3f' % boxLength)
            if nCount%3 == 0:
                outFile.write('  90.00  90.00  90.00 P  1  1  1    1\n')
    if len(words) == 9: # lipids+water system
      n=int(words[0]) + nExtra # atom identifier
      t=int(words[1]) # type identifier
      m=int(words[2]) # molecule identifier
      x=float(words[3])
      y=float(words[4])
      z=float(words[5])
      mux=float(words[6])
      muy=float(words[7])
      muz=float(words[8])
      muMag=sqrt(mux*mux+muy*muy+muz*muz)
      if muMag>0.0:
        watType = t
        # compute coordinate of dipole's "+" tip:
        xPlus=x+SCALE*mux/muMag
        yPlus=y+SCALE*muy/muMag
        zPlus=z+SCALE*muz/muMag
        # compute coordinate of dipole's "-" tail:
        xMinus=x-SCALE*mux/muMag
        yMinus=y-SCALE*muy/muMag
        zMinus=z-SCALE*muz/muMag
        outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                      %(n,"NN NEG",m,xMinus,yMinus,zMinus ))
        outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                      %(n+1,"NP POS",m,xPlus,yPlus,zPlus ))
        nExtra=nExtra+1
      else:
        if t>watType: # if this is an ion
          outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                        %(n,"I     ",m,xMinus,yMinus,zMinus ))        
print "File \'elbawater.pdb\' successfully generated."
