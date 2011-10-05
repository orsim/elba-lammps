#!/usr/bin/env python

# Script:  trj2pdb.py
# Purpose: this script reads a LAMMPS ".trj" trajectory (dump) file and 
#          converts it into a trajectory file "trajectory.pdb", which can
#          be read in VMD.
#          Point-dipoles in LAMMPS ".trj" files are described by mass center
#          coordinates and orientations (x-, y-, and z- projections of the
#          dipole vectors). Since VMD does not allow visualization of point-
#          dipoles, this script uses the orientation information to convert
#          every dipolar atom into two atoms representing the "+" and "-"
#          ends of the original dipole. Such two atoms are spaced 1 Angstrom
#          apart. 
# Syntax:  trj2pdb.py dump.trj

import sys, string, linecache
from math import sqrt

inFileName = sys.argv[1]
inFile = open(inFileName, "r")

print "Processing file %s ..." % inFileName
line = linecache.getline(inFileName, 4)
words = string.split(line)
nAtoms = int(words[0])
print "Number of atoms: %d" % nAtoms

outFile = open("trajectory.pdb", "w")
lines = inFile.readlines()

HALF_ANGSTROM=0.5; # [Angstrom] scaling factor
nLine=0; # line counter
nCount=0; # counter
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
            outFile.write('CRYST1')
        else:
            nCount=nCount+1
            negCoord=float(words[0])
            plusCoord=float(words[1])
            boxLength=plusCoord-negCoord
            outFile.write('%9.3f' % boxLength)
            if nCount%3 == 0:
                outFile.write('  90.00  90.00  90.00 P  1  1  1    1\n')
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
        xPlus=x+HALF_ANGSTROM*mux/muMag
        yPlus=y+HALF_ANGSTROM*muy/muMag
        zPlus=z+HALF_ANGSTROM*muz/muMag
        # compute coordinate of dipole's "-" tail:
        xMinus=x-HALF_ANGSTROM*mux/muMag
        yMinus=y-HALF_ANGSTROM*muy/muMag
        zMinus=z-HALF_ANGSTROM*muz/muMag
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
            xPlus=x+HALF_ANGSTROM*mux/muMag
            yPlus=y+HALF_ANGSTROM*muy/muMag
            zPlus=z+HALF_ANGSTROM*muz/muMag
            # compute coordinate of dipole's "-" tail:
            xMinus=x-HALF_ANGSTROM*mux/muMag
            yMinus=y-HALF_ANGSTROM*muy/muMag
            zMinus=z-HALF_ANGSTROM*muz/muMag
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
