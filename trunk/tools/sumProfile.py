#!/usr/bin/env python

# Script: sumProfile.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Reads two or more profiles and calculates sum
# Syntax: sumProfile.py fileName1 fileName2 ... fileNameN
# Example: sumProfile.py eppHead.dat eppWat.dat > eppSum.dat

import sys,string

if len(sys.argv) <= 2:
  print "Syntax: sumProfile.py fileName1 fileName2 ... fileNameN"
  sys.exit()

nInputProfiles = len(sys.argv[1:])
coordinates = []
values = []

# store all data from input profiles:
for i in range(0, nInputProfiles):
    inFile = open(sys.argv[i+1],'r')
    lines = inFile.readlines()
    nValues_singleProf = len(lines)
    for line in lines:
        words = string.split(line)
        coordinates.append(float(words[0]))    
        values.append(float(words[1]))

# evaluate profile sum:
for i in range(0, nValues_singleProf):
    coord=coordinates[i]
    sumVal=0
    for j in range(0, nInputProfiles):
        sumVal = sumVal + values[i+j*nValues_singleProf]
    print coord, sumVal
    
