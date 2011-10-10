#!/usr/bin/env python

# Script: lpp2st.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Integrates a lateral pressure profile to get the corresponding
#          surface tension per monolayer
# Syntax: lpp2st.py inputFile
#         inputFile = coordinate | lpp(z)
# Notes: - coordinates (first col of input file) is in Angstrom
#        - lateral pressure Pt-Pn (second col of input file) is in Atm 
# Example: lpp2st.py lpp.dat 
# References: - Alejandre et al, J Chem Phys 102, 4574 (1995)

import sys,os,string, linecache

def Integral( iStart, iEnd, f, delta ):
  fInt = 0.0
  for i in range( iStart, iEnd ): # i = iStart, ..., iEnd-1
    fInt += f[i]
  return fInt*delta

cFac = 1.01325e-2 # conversion factor: atm*A = 1.01325e-5 mN/m

if len(sys.argv) != 2:
  print "Syntax: lpp2st.py inputFile"
  sys.exit()

inFileName = sys.argv[1]
inFile = open(inFileName, "r")
lines = inFile.readlines()
inFile.close()

# find slab thickness (delta), #slabs and box size:
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

st = - 0.5 * Integral( 0, nSlabs, lp, delta ) * cFac
print 'Surface tension per monolayer from entire lpp: %f mN/m' % st

st = - Integral( 0, nSlabs/2, lp, delta ) * cFac
print 'Surface tension on \'-\' monolayer: %f mN/m' % st

st = - Integral( nSlabs/2, nSlabs, lp, delta ) * cFac
print 'Surface tension on \'+\' monolayer: %f mN/m' % st


