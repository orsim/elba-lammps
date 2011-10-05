#!/usr/bin/env python

# Script:  box.py
# Purpose: this script reads a LAMMPS ".trj" trajectory (dump) file and 
#          extracts information on the simulation box.
# Output: 'vol.dat': box volume [angstrom^3]
#         'numDens.dat': number density [1/angstrom^3] 
#         'watDens.dat': water density [g/cm^3]
#         'xyArea.dat': area of xy plane [angstrom^2]
#         'xzArea.dat': area of xz plane [angstrom^2]
#         'yzArea.dat': area of yz plane [angstrom^2]
# Usage:   python box.py dump.trj

import sys, string, linecache
from math import sqrt

A3_IN_cm3 = 1e-24; # conversion factor: Angstrom^3 -> cm^3
#watMass_g = 3e-23; # water molecular mass in grams
watMass_g = 18.01528/6.0221367e23; # water molecular mass in grams

inFileName = sys.argv[1]
inFile = open(inFileName, "r")

print "Processing file %s ..." % inFileName
line = linecache.getline(inFileName, 4)
words = string.split(line)
nAtoms = int(words[0])
print "Number of atoms: %d" % nAtoms

outVolFile = open("vol.dat", "w")
outXYFile = open("xyArea.dat", "w")
outXZFile = open("xzArea.dat", "w")
outYZFile = open("yzArea.dat", "w")
outNumDensFile = open("numDens.dat", "w")
outWatDensFile = open("watDens.dat", "w")

lineNum=0; # line counter
snapNum=1;
boxLine=5; # first instance of 'ITEM: BOX BOUNDS pp pp pp'
lines = inFile.readlines()
for line in lines:
    lineNum = lineNum + 1
    words = string.split(line)
    if len(words) == 6:
        if words[1]=="BOX":
            boxLine=lineNum
    if lineNum==boxLine+1:
        xSize = float(words[1]) - float(words[0])
    if lineNum==boxLine+2:
        ySize = float(words[1]) - float(words[0])
    if lineNum==boxLine+3:
        zSize = float(words[1]) - float(words[0])
        volume = xSize * ySize * zSize
        xyArea = xSize * ySize
        xzArea = xSize * zSize
        yzArea = ySize * zSize
        outVolFile.write("%f\n" % volume)
        outXYFile.write("%f\n" % xyArea)
        outXZFile.write("%f\n" % xzArea)
        outYZFile.write("%f\n" % yzArea)
        outNumDensFile.write("%f\n" % (nAtoms/volume))
        outWatDensFile.write("%f\n" % (nAtoms*watMass_g/(volume*A3_IN_cm3)))
        
