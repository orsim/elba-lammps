#!/usr/bin/env python

# Script: calcHeatOfVap.py
# Purpose: Calculates (molar) heat of vaporization from intermolecular
#          potential energy per mole U and temperature T
# Syntax: calcHeatOfVap.py U T
#         U = intermolecular potential energy per mole [kcal/mol]
#         T = temperature [K]
# Examples: calcHeatOfVap.py -10 298
# Notes: Benchmark values for the heat of vaporization at 298 K:
#        i) 11.0 kcal/mol from 'classical calculation' [1,2]
#        ii) 10.52 kcal/mol from experiment [3]
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# References: [1] Guillot, J Mol Liq 101, 219 (2002)
#             [2] Guillot & Guissani, J Chem Phys 114, 6720 (2001)
#             [3] Haar, Gallagher & Kell, NBS/NRC Steam tables,
#                 Hemisphere, Washington DC, 1984

import sys,os,string

if len(sys.argv) != 3:
  print "Syntax: calcHeatOfVap.py U T"
  sys.exit()

U = float(sys.argv[1])
T = float(sys.argv[2])

R = 1.987e-3 # gas constant in kcal/(mol K)

print "%f kcal/mol" % ( -U + R*T )
