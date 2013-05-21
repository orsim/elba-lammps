#!/usr/bin/env python

# Script: lpp2gamma.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Integrates a lateral pressure profile to get the corresponding
#          surface tension (gamma) per monolayer
# Syntax: lpp2gamma.py inputFile
# Example: lpp2gamma.py lpp.dat 
# Notes: - inputFile format = coordinate | lpp(z)
#        - Coordinates (first col of input file) in Angstrom
#        - Lateral pressure Pt-Pn (second col of input file) in Atm 
#        - Coordinates' origin = bilayer center
#        - Units: mN/m = dyn/cm
# Reference: Alejandre et al, J Chem Phys 102, 4574 (1995)

import sys,os,string

if len(sys.argv) != 2:
  print "Syntax: lpp2gamma.py inputFile"
  sys.exit()

def Integral( iStart, iEnd, f, delta ):
  fInt = 0.0
  for i in range( iStart, iEnd ): # i = iStart, ..., iEnd-1
    fInt += f[i]
  return fInt*delta

atmA__in__mN_m = 0.01013 # conversion factor: atm*A = 0.01013 mN/m

inFileName = sys.argv[1]
inFile = open(inFileName, "r")
lines = inFile.readlines()
inFile.close()

# find slab thickness (delta), number of slabs, and box size:
coord = [] # coordinate
lp = [] # lateral pressure
nSlabs = 0
for line in lines:
  words = string.split(line)
  if len(words) == 2:
    nSlabs = nSlabs + 1
    coord.append(float(words[0]))
    lp.append(float(words[1]))
zBox = coord[nSlabs-1] - coord[0]
delta = zBox / ( nSlabs - 1 )

st = - Integral( 0, nSlabs/2, lp, delta ) * atmA__in__mN_m
print 'Surface tension on \'-\' monolayer: %f mN/m' % st

st = - Integral( nSlabs/2, nSlabs, lp, delta ) * atmA__in__mN_m
print 'Surface tension on \'+\' monolayer: %f mN/m' % st

st = - 0.5 * Integral( 0, nSlabs, lp, delta ) * atmA__in__mN_m
print 'Avg surf tens per monolayer from entire lpp: %f mN/m' % st


