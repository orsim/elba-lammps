#!/usr/bin/env python

# Script: rdf2dat.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Reads a LAMMPS output file containing radial distribution
#          function g(r) data and sorts such data in |r|g| format
# Syntax: rdf2dat.py inputFile 
# Examples: - rdf2dat.py wat.rdf > rdf.dat

import sys, string

if len(sys.argv) != 2:
  print "Syntax: rdf2dat.py inputFile"
  sys.exit()

inFileName = sys.argv[1]
inFile = open(inFileName, "r")
lines = inFile.readlines()
inFile.close()

# read input data:
for line in lines:
    if line[0] != '#': # ignore comments
        words = string.split(line)
        if len(words) != 2: # ignore "timestep/rows" values
          print words[1], words[2] # print coord, rdf
