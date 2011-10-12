#!/usr/bin/env python

# Script: lpp2curv.py
# Purpose: Process a lateral pressure profile to extract curvature
#          elastic parameters (per monolayer)
# Syntax: lpp2curv.py inputFile
# Example: lpp2curv.py lpp.dat 
# Notes: - inputFile format = coordinate | lpp(z)
#        - Coordinates (first col of input file) in Angstrom
#        - Lateral pressure Pt-Pn (second col of input file) in Atm
#        - Coordinates' origin = bilayer center
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
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

tau1m_z = [] # fist int mom profile
tau2m_z = [] # second int mom profile
outFile1 = open('tau1m.dat', 'w')
outFile2 = open('tau2m.dat', 'w')

# calc 1st integral moment of lpp across 'negative' monolayer:
tau1mNeg = 0 # initialize 1st mom integral sum
tau2mNeg = 0 # initialize 2nd mom integral sum
for i in range( nSlabs/2 - 1, -1, -1 ): # i = nSlabs/2 - 1, ..., 0
  tau1mNeg += (-coord[i]) * lp[i] * delta
  tau2mNeg += coord[i]**2 * lp[i] * delta
  tau1m_z.insert(0,tau1mNeg) # insert at the front
  tau2m_z.insert(0,tau2mNeg) # insert at the front
tau1mNeg = tau1mNeg * atmA2__in__J_nm
tau2mNeg = tau2mNeg * atmA3__in__J
print ( '1st int mom for \'-\' monolayer: tau1m = %e J/nm = %f kBT/nm'
        % ( tau1mNeg, tau1mNeg * J__in__kBTroom ) )
print ( '2nd int mom for \'-\' monolayer: tau2m = %e J = %f kBT'
        % ( tau2mNeg, tau2mNeg * J__in__kBTroom ) )

# sort and dump 'negative' profiles:
for i in range( 0, nSlabs/2 ):
  outFile1.write( '%f %f\n' % (coord[i],-tau1m_z[i]) )
  outFile2.write( '%f %f\n' % (coord[i], tau2m_z[i]) )

# calc 1st integral moment of lpp across 'positive' monolayer:
tau1mPos = 0 # initialize 1st mom integral sum
tau2mPos = 0 # initialize 2nd mom integral sum
for i in range( nSlabs/2, nSlabs ): # i = nSlabs/2, ..., nSlabs-1
  tau1mPos += coord[i] * lp[i] * delta
  tau2mPos += coord[i]**2 * lp[i] * delta
  outFile1.write( '%f %f\n' % (coord[i],tau1mPos) )
  outFile2.write( '%f %f\n' % (coord[i],tau2mPos) )
tau1mPos = tau1mPos * atmA2__in__J_nm
tau2mPos = tau2mPos * atmA3__in__J
print ( '1st int mom for \'+\' monolayer: tau1m = %e J/nm = %f kBT/nm'
        % ( tau1mPos, tau1mPos * J__in__kBTroom ) )
print ( '2nd int mom for \'+\' monolayer: tau2m = %e J = %f kBT'
        % ( tau2mPos, tau2mPos * J__in__kBTroom ) )

# calc & print averages:
tau1mAvg = 0.5 * ( tau1mNeg + tau1mPos )
print ( 'tau1m averaged over the two monolayers = %e J/nm = %f kBT/nm' 
        % ( tau1mAvg, tau1mAvg * J__in__kBTroom ) )
tau2mAvg = 0.5 * ( tau2mNeg + tau2mPos )
print ( 'tau2m averaged over the two monolayers = %e J = %f kBT'
        % ( tau2mAvg, tau2mAvg * J__in__kBTroom ) )
