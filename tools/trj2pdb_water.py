## Copyright (C) 2011 Mario Orsi
## This file is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at your option)
## any later version. This file is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
## Public License for more details: <http://www.gnu.org/licenses/>. 

#############################################################################
# This script reads the LAMMPS dump file "dump.lammpstrj" and outputs the
# trajectory file "trajectory.pdb" which can be read in VMD.
#
# Dipoles in "dump.lammpstrj" are described by mass center coordinates and
# orientations. Since VMD does not (currently) allow visualization of dipoles,
# this script uses the orientation information to add an extra site to each
# atom; such an extra site corresponds to the "+" end of the dipole
#
# USAGE: python lammpstrj2pdb.py

import sys, string, linecache
from math import sqrt

print "Processing file dump.lammpstrj ..."
line = linecache.getline('dump.lammpstrj', 4)
words = string.split(line)
nSites = int(words[0])
print "Number of sites: %d" % nSites

outFile = open("trajectory.pdb", "w")

inFile = open("dump.lammpstrj", "r")
lines = inFile.readlines()

print "Setting dipole scaling factor assuming real units..."
s=2; # scaling factor
#s=s/30; # rescale for reduced units - comment this line for real units
nLine=0; # line counter

print "Start scanning coordinates and orientations..."
for line in lines:
    nLine=nLine+1
    words = string.split(line)
    if len(words) == 2:
        if words[1]=="TIMESTEP":
            if nLine!=1:
                outFile.write('ENDMDL\n')
            outFile.write('MODEL\n')
    if len(words) == 8:
        n=int(words[0]) # atom identifier
        t=int(words[1]) # type identifier
        x=float(words[2])
        y=float(words[3])
        z=float(words[4])
        mux=float(words[5])
        muy=float(words[6])
        muz=float(words[7])
        xh=x+s*mux
        yh=y+s*muy
        zh=z+s*muz
        outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                      %(n,"OW WAT",n,x,y,z ))
        nPlus=n+nSites # index of added hydrogen ('+' end of dipole)
        outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                      %(nPlus,"HW WAT",n,xh,yh,zh ))
        
print "Done scanning - end of script."
