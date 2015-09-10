#!/usr/bin/env python

# Script: calcVolComp.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Reads a file containing the time evolution of the simulation
#          box volume and calculates the corresponding volume
#          compressibility (bulk) modulus
# Syntax: python calcVolComp.py inputFile temperature
# Notes: - input file is assumed to have 2-column "step | volume" format
#        - input volume is in Angstrom^3
#        - input file 1st column (step) is not used
# Example: python calcVolComp.py boxVol.dat 303.15
# Reference: - Orsi et al, J Phys Condens Matter 22: 155106 (2010),
#              section 5.1.2

import sys, string
from math import sqrt

kB = 1.3806505e-23 # Boltzmann constant [J / K]
A3_in_nm3 = 1e-3 # 1 Angstrom^3 = 10^(-3) nm^3
A6_in_nm6 = 1e-6 # 1 Angstrom^6 = 10^(-6) nm^6
J_nm3__in__kbar = 1e19 # 1 J/nm^3 = 10^19 kbar

if len(sys.argv) != 3:
  print "Syntax: python calcVolComp.py inputFile temperature"
  sys.exit()

inFileName = sys.argv[1]
T = float(sys.argv[2])
print "\nT = %.2f K\n" % T
inFile = open(inFileName, "r")
lines = inFile.readlines()
inFile.close()

dataCounter = 0

volSum = squareVolSum = 0
volSum1 = squareVolSum1 = 0
volSum2 = squareVolSum2 = 0

nData = nData1 = nData2 = 0 # counter for number of data (measurements)

for line in lines:
    if line[0] != '#': # ignore comments
        nData += 1

for line in lines: 
    if line[0] != '#': # ignore comments
        dataCounter = dataCounter + 1
        words = string.split( line )
        vol = string.atof(words[1]);
        volSum += vol
        squareVolSum += vol**2
        if dataCounter <= nData/2 : # first half of data
            nData1 += 1
            volSum1 += vol
            squareVolSum1 += vol**2
        else: # second half of data
            nData2 += 1
            volSum2 += vol
            squareVolSum2 += vol**2
            
# check:
if nData != nData1+nData2:
    print "ERROR: wrong scanning of data - check script"

# Calc modulus considering all data:
meanVol = volSum / nData;
print "Mean box volume = %4.2f A^3" % meanVol
meanSquaredVolFluct = squareVolSum / nData - meanVol**2
#print "meanSquaredVolFluct = %6.3f A^6" % meanSquaredVolFluct
# convert A -> nm:
meanVol *= A3_in_nm3 
meanSquaredVolFluct *= A6_in_nm6
# computing modulus:
KV = kB*T * meanVol / meanSquaredVolFluct # [ J / nm^2 ]
# convert J/nm^2 -> kbar:
KV *= J_nm3__in__kbar
print "KV_total = %.2f kbar\n" % KV

# Calc modulus considering only first half of data:
meanVol1 = volSum1 / nData1;
#print "meanVol1 = %4.2f A^3" % meanVol1
#print "squareVolSum1 / nData1= %4.2f A^3" % (squareVolSum1/ nData1)
#print "meanVol1**2 = %4.2f A^3" % (meanVol1**2)
meanSquaredVolFluct1 = squareVolSum1 / nData1 - meanVol1**2
#print "meanSquaredVolFluct1 = %6.3f A^6" % meanSquaredVolFluct1
# convert A -> nm:
meanVol1 *= A3_in_nm3 
meanSquaredVolFluct1 *= A6_in_nm6
# computing modulus:
KV1 = kB*T * meanVol1 / meanSquaredVolFluct1 # [ J / nm^2 ]
# conversion considering that J/nm^2 = 10^21 kbar
KV1 *= J_nm3__in__kbar
print "KV_1 = %.2f kbar\n" % ( KV1 )

# Calc modulus considering only second half of data:
meanVol2 = volSum2 / nData2;
#print "meanVol2 = %4.2f A^3" % meanVol2
#print "squareVolSum2/ nData2 = %4.2f A^3" % (squareVolSum2/ nData2)
#print "meanVol2**2 = %4.2f A^3" % (meanVol2**2)
meanSquaredVolFluct2 = squareVolSum2 / nData2 - meanVol2**2
#print "meanSquaredVolFluct2 = %6.3f A^6" % meanSquaredVolFluct2
# convert A -> nm:
meanVol2 *= A3_in_nm3 
meanSquaredVolFluct2 *= A6_in_nm6
# computing modulus:
KV2 = kB*T * meanVol2 / meanSquaredVolFluct2 # [ J / nm^2 ]
# conversion considering that J/nm^2 = 10^21 kbar
KV2 *= J_nm3__in__kbar
print "KV_2 = %.2f kbar\n" % ( KV2 )

# Calc average:
KV_avg12 = 0.5 * ( KV1 + KV2 )
print "KV_avg12 = %.2f kbar" % KV_avg12
standardDeviation = sqrt( ( KV1 - KV_avg12 )**2 + ( KV2 - KV_avg12 )**2 )
standardError = standardDeviation / sqrt(2)
print "Standard error = %.2f kbar\n" % standardError
