#!/usr/bin/env python

# Script: calcHeatCap.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Calculates (molar) heat capacity at constant pressure
#          from total energy per molecule E
# Syntax: calcHeatCap.py E1 E2 T1 T2
#         E1 = total energy per molecule [kcal/mol] at T1
#         E2 = total energy per molecule [kcal/mol] at T2
#         T1 = temperature [K]
#         T2 = temperature [K]
# Examples: calcHeatCap.py -10 -12 298 318
# Notes: - The calculation includes quantum corrections [1] 
#        - Benchmark value for the heat capacity is ~18 kcal/mol [3]
# References: [1] Glattli et al, J Chem Phys 2002, 116, 9811
#             [3] Weast, Handbook of chemistry and physics, 61st ed, CRC
#                 Boca Raton, FL, 1980

import sys,os

if len(sys.argv) != 5:
  print "Syntax: calcHeatCap.py E1 E2 T1 T2"
  sys.exit()

E1 = float(sys.argv[1])
E2 = float(sys.argv[2])
T1 = float(sys.argv[3])
T2 = float(sys.argv[4])

Q = -2.22 # [kcal/mol] quantum contribution (see ref [1] above)

print "%f kcal/mol" % ( (E1-E2)/(T1-T2) + Q )
