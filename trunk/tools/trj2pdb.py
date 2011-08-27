# Copyright (C) 2011 Mario Orsi
# This file is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. This file is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details: <http://www.gnu.org/licenses/>. 

#---------------------------------------------------------------------------
# This script reads the LAMMPS dump file "dump.trj" and outputs the
# trajectory file "trajectory.pdb", which can be read in VMD.
#
# Dipoles in "dump.trj" are described by mass center coordinates and
# orientations (x-, y-, and z- projections of the dipole vectors). Since
# VMD does not (currently) allow visualization of dipoles, this script uses
# the orientation information to convert every dipolar atom into two atoms
# representing the "+" and "-" ends of the original dipole. Such two atoms
# are spaced 1 Angstrom apart. 
#
# USAGE: $python trj2pdb.py
#---------------------------------------------------------------------------

import sys, string, linecache
from math import sqrt

inFile = open("dump.trj", "r")

print "Processing file dump.trj ..."
line = linecache.getline('dump.trj', 4)
words = string.split(line)
nSites = int(words[0])
print "Number of sites: %d" % nSites

outFile = open("trajectory.pdb", "w")

lines = inFile.readlines()

s=0.5; # [Angstrom] scaling factor
nLine=0; # line counter
nExtra=0; # extra atoms used to represent dipoles

print "Start scanning coordinates and orientations..."
for line in lines:
    nLine=nLine+1
    words = string.split(line)
    if len(words) == 2:
        if words[1]=="TIMESTEP":
            if nLine!=1:
                outFile.write('ENDMDL\n')
                nExtra=0;
            outFile.write('MODEL\n')
    if len(words) == 8: # water-only system
        n=int(words[0]) + nExtra # atom identifier
        t=int(words[1]) # type identifier
        x=float(words[2])
        y=float(words[3])
        z=float(words[4])
        mux=float(words[5])
        muy=float(words[6])
        muz=float(words[7])
        muMag=sqrt(mux*mux+muy*muy+muz*muz)
        # compute coordinate of dipole's "+" tip:
        xPlus=x+s*mux/muMag
        yPlus=y+s*muy/muMag
        zPlus=z+s*muz/muMag
        # compute coordinate of dipole's "-" tail:
        xMinus=x-s*mux/muMag
        yMinus=y-s*muy/muMag
        zMinus=z-s*muz/muMag
        outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                      %(n,"OW WAT",nExtra+1,xMinus,yMinus,zMinus ))
        outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                      %(n+1,"HW WAT",nExtra+1,xPlus,yPlus,zPlus ))
        nExtra=nExtra+1
    if len(words) == 9: # lipids+water system
        n=int(words[0]) + nExtra # atom identifier
        t=int(words[1]) # type identifier
        m=float(words[2]) # molecule identifier
        x=float(words[3])
        y=float(words[4])
        z=float(words[5])
        mux=float(words[6])
        muy=float(words[7])
        muz=float(words[8])
        muMag=sqrt(mux*mux+muy*muy+muz*muz)
        if muMag>0.0:
            # compute coordinate of dipole's "+" tip:
            xPlus=x+s*mux/muMag
            yPlus=y+s*muy/muMag
            zPlus=z+s*muz/muMag
            # compute coordinate of dipole's "-" tail:
            xMinus=x-s*mux/muMag
            yMinus=y-s*muy/muMag
            zMinus=z-s*muz/muMag
        if t == 1: # water type
            outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                          %(n,"OW WAT",m,xMinus,yMinus,zMinus ))
            outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                          %(n+1,"HW WAT",m,xPlus,yPlus,zPlus ))
            nExtra=nExtra+1
        elif t == 2: # choline type
            outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                          %(n,"CHO LIP",m,x,y,z ))
        elif t == 3: # phosphate type
            outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                          %(n,"PHO LIP",m,x,y,z ))
        elif t == 4 or t == 5: # glycerol or ester type
            if t == 4: # glycerol type
                outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                              %(n,"Gp LIP",m,xPlus,yPlus,zPlus ))
                outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                              %(n+1,"Gm LIP",m,xMinus,yMinus,zMinus ))
                nExtra=nExtra+1
            if t == 5: # ester type
                outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                              %(n,"Em LIP",m,xMinus,yMinus,zMinus ))
                outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                              %(n+1,"Ep LIP",m,xPlus,yPlus,zPlus ))
                nExtra=nExtra+1
        elif t == 6: # tail type
            outFile.write('ATOM%7d%9s%6d%12.3f%8.3f%8.3f\n'
                          %(n,"CH2 LIP",m,x,y,z ))
        
print "Done scanning - end of script."
