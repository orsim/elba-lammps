#!/usr/bin/env python

# Script: lpp2st.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Integrates a lateral pressure profile to get the corresponding
#          surface tension per monolayer
# Syntax: lpp2st.py inputFile
#         inputFile = coordinate | lpp(z)
# Notes: - coordinates (first col of input file) is in Angstrom
#        - lateral pressure (second col of input file) is in Atm 
# Example: lpp2st.py lpp.dat 
# References: - Harries & Ben-Shaul, J Chem Phys 106, 1609 (1997)

import sys,os,string, linecache

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
    coord.append(words[0])
    lp.append(words[1])
zBox = float(coord[nSlabs-1]) - float(coord[0])
delta = zBox / ( nSlabs - 1 )

# integrate lpp across whole box
lpInt = 0.5 * ( float(lp[0]) + float(lp[nSlabs-1]) )
for i in range( 1, nSlabs-1 ): # i = 1, ..., nSlabs-2
  lpInt += float( lp[i] )
st = - 0.5 * lpInt * delta * cFac
print 'Surface tension per monolayer from entire lpp: %f mN/m' % st

nSlabsHalf = nSlabs / 2;

# integrate lpp across first half-box
lpInt = 0.5 * ( float(lp[0]) + float(lp[nSlabsHalf-1]) )
for i in range( 1, nSlabsHalf-1 ): # i = 1, ..., nSlabsHalf-2
  lpInt += float( lp[i] )
st = - lpInt * delta * cFac
print 'Surface tension on \'-\' monolayer: %f mN/m' % st

# integrate lpp across second half-box
lpInt = 0.5 * ( float(lp[nSlabsHalf]) + float(lp[nSlabs-1]) )
for i in range( nSlabsHalf+1, nSlabs-1 ): # i = nSlabsHalf+1, ..., nSlabs-2
  lpInt += float( lp[i] )
st = - lpInt * delta * cFac
print 'Surface tension on \'+\' monolayer: %f mN/m' % st


