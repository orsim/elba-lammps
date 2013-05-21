#!/usr/bin/env python

# Script: numDens2elecDens.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Reads a file containing a number density distribution
#          and calculates the corresponding electron density
# Syntax: numDens2elecDens.py inputFile nElectrons
# Example: numDens2elecDens.py numDensChol.zProfile 50 > edpChol.dat
# Notes:  - inputFile = LAMMPS output file generated by fix ave/spatial
#         - nElectrons = equivalent number of electrons per CG site
# References: - Orsi & Essex, Faraday Discuss 161, 249 (2013)
#             - Orsi & Essex, PLoS ONE 6, e28637 (2011)
#             - Orsi et al, J Phys Condens Matter 22 (2010) 155106,
#             section 5.1.4

import sys,os,string

if len(sys.argv) != 3:
  print "Syntax: numDens2elecDens.py inputFile nElectrons"
  sys.exit()

inFileName = sys.argv[1]
nElectrons = float(sys.argv[2])
inFile = open(inFileName, "r")
lines = inFile.readlines()
inFile.close()

# calculate and electron density:
for line in lines:
    if line[0] != '#': # ignore comments
        words = string.split(line)
        if len(words) == 4:
            coord = float(words[1])
            electronDensity = float(words[3]) * nElectrons
            print coord, electronDensity # [A, e/A^3]