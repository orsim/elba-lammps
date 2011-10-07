#!/usr/bin/env python

# Script: lpp2curv.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Process a lateral pressure profile to extract curvature
#          elastic parameters (per monolayer)
# Syntax: lpp2curv.py inputFile
#         inputFile = coordinate | lpp(z)
# Notes: - coordinates (first col of input file) in Angstrom
#        - lateral pressure Pt-Pn (second col of input file) in Atm 
# Example: lpp2curv.py lpp.dat 
# References: - Orsi et al, J Phys Condens Matter 22, 155106 (2010),
#               section 5.3
#             - Cantor, Biophys J 80, 2284 (2001)


import sys,os,string, linecache

atmA2__in__J_nm = 1.01325e-24 # atm*A^2 = 1.01325e-24 J/nm
J__in__kBTroom = 1.0/4.0453e-21 # kB*Troom = 4.0453e-21 J (Troom = 293 K)

cFac = 0.00025048

if len(sys.argv) != 2:
  print "Syntax: lpp2curv.py inputFile"
  sys.exit()

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

nSlabsHalf = nSlabs / 2;

tau1m_z = [] # fist int mom profile 
outFile = open('tau1m.dat', 'w')

# calc 1st integral moment of lpp across 'negative' monolayer:
z0 = -coord[nSlabsHalf-1] # coord at bilayer center
zh = -coord[0] # coord at box edge (int is over half-box)
tau1m_0 = 0.5 * delta * z0*lp[nSlabsHalf-1]
tau1m_h = 0.5 * delta * zh*lp[0]
tau1m = tau1m_0 + tau1m_h 
tau1m_z.append( tau1m_0 )
outFile.write( '%f %f\n' % ( -zh, tau1m_h ) )
for i in range( 1, nSlabsHalf-1 ): # i = 1, ..., nSlabsHalf-2
  tau1m += - delta * coord[i] * lp[i]
  tau1m_z.append(tau1m)
  outFile.write( '%f %f\n' % (coord[i],tau1m) )
tau1m = tau1m * atmA2__in__J_nm
outFile.write('%f %f\n' % (-z0, tau1m_0))
print ( '1st int mom for \'-\' monolayer: tau1m = %e J/nm = %f kBT/nm'
        % ( tau1m, tau1m * J__in__kBTroom ) )

# calc 1st integral moment of lpp across 'positive' monolayer:
z0 = coord[nSlabsHalf] # coord at bilayer center
zh = coord[nSlabs-1] # coord at box edge (int is over half-box)
tau1m_0 = 0.5 * delta * z0*lp[nSlabsHalf]
tau1m_h = 0.5 * delta * zh*lp[nSlabs-1]
tau1m = tau1m_0 + tau1m_h 
tau1m_z.append( tau1m_0 )
outFile.write( '%f %f\n' % ( z0, tau1m_0) )
for i in range( nSlabsHalf+1, nSlabs-1 ): # i = nSlabsHalf+1, ..., nSlabs-2
  tau1m += delta * coord[i] * lp[i] 
  tau1m_z.append(tau1m)
  outFile.write( '%f %f\n' % (coord[i],tau1m) )
tau1m = tau1m * delta * atmA2__in__J_nm
outFile.write( '%f %f\n' % (zh, tau1m_h) )
print ( '1st int mom for \'+\' monolayer: tau1m = %e J/nm = %f kBT/nm'
        % ( tau1m, tau1m * J__in__kBTroom ) )


