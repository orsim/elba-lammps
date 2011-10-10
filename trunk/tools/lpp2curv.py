#!/usr/bin/env python

# Script: lpp2curv.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Process a lateral pressure profile to extract curvature
#          elastic parameters (per monolayer)
# Syntax: lpp2curv.py inputFile
#         inputFile = coordinate | lpp(z)
# Notes: - Coordinates (first col of input file) in Angstrom
#        - Lateral pressure Pt-Pn (second col of input file) in Atm
#        - Coordinates' origin = bilayer center
#        - Not tested on odd numbers of slabs
# Example: lpp2curv.py lpp.dat 
# References: - Orsi & Essex, PLoS ONE, submitted
#             - Orsi et al, J Phys: Condens Matter 22, 155106 (2010)
#             - Cantor, Biophys J 80, 2284 (2001)
#             - Ben-Shaul, in 'Structure and dynamics of membranes',
#               eds Lipowsky & Sackmann, Elsevier (1995)  

import sys,os,string

if len(sys.argv) != 2:
  print "Syntax: lpp2curv.py inputFile"
  sys.exit()

atmA2__in__J_nm = 1.01325e-24 # atm*A^2 = 1.01325e-24 J/nm
atmA3__in__J = 1.01325e-25 # atm*A^3 = 1.01325e-25 J
J__in__kBTroom = 1.0/4.0453e-21 # kB*Troom = 4.0453e-21 J (Troom=293 K)

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
tau2m_z = [] # second int mom profile
outFile1 = open('tau1m.dat', 'w')
outFile2 = open('tau2m.dat', 'w')

# calc 1st integral moment of lpp across 'negative' monolayer:
tau1m = 0 # initialize 1st mom integral sum
tau2m = 0 # initialize 2nd mom integral sum
for i in range( nSlabsHalf-1, -1, -1 ): # i = nSlabsHalf-1, ..., 0
  tau1m += (-coord[i]) * lp[i] * delta
  tau2m += coord[i]**2 * lp[i] * delta
  tau1m_z.insert(0,tau1m) # insert at the front
  tau2m_z.insert(0,tau2m) # insert at the front
tau1m = tau1m * atmA2__in__J_nm
tau2m = tau2m * atmA3__in__J
print ( '1st int mom for \'-\' monolayer: tau1m = %e J/nm = %f kBT/nm'
        % ( tau1m, tau1m * J__in__kBTroom ) )
print ( '2nd int mom for \'-\' monolayer: tau2m = %e J = %f kBT'
        % ( tau2m, tau2m * J__in__kBTroom ) )

# sort and dump 'negative' profiles:
for i in range( 0, nSlabsHalf ):
  outFile1.write( '%f %f\n' % (coord[i],-tau1m_z[i]) )
  outFile2.write( '%f %f\n' % (coord[i], tau2m_z[i]) )

# calc 1st integral moment of lpp across 'positive' monolayer:
tau1m = 0 # initialize 1st mom integral sum
tau2m = 0 # initialize 2nd mom integral sum
for i in range( nSlabsHalf, nSlabs ): # i = nSlabsHalf, ..., nSlabs-1
  tau1m += coord[i] * lp[i] * delta
  tau2m += coord[i]**2 * lp[i] * delta
  outFile1.write( '%f %f\n' % (coord[i],tau1m) )
  outFile2.write( '%f %f\n' % (coord[i],tau2m) )
tau1m = tau1m * atmA2__in__J_nm
tau2m = tau2m * atmA3__in__J
print ( '1st int mom for \'+\' monolayer: tau1m = %e J/nm = %f kBT/nm'
        % ( tau1m, tau1m * J__in__kBTroom ) )
print ( '2nd int mom for \'+\' monolayer: tau2m = %e J = %f kBT'
        % ( tau2m, tau2m * J__in__kBTroom ) )


