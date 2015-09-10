#!/usr/bin/env python

# Script: calcHeatCap.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Calculates (molar) heat capacity at constant pressure
#          from total energy per molecule E
# Syntax: python calcHeatCap.py E1 E2 T1 T2
#         E1 = total energy per molecule [kcal/mol] at T1
#         E2 = total energy per molecule [kcal/mol] at T2
#         T1 = temperature [K]
#         T2 = temperature [K]
# Example: python calcHeatCap.py -10 -12 298 318
# Notes: - The calculation includes quantum corrections [1] 
#        - Benchmark value for the heat capacity is ~18 kcal/mol [3]
# References: [1] Glattli et al, J Chem Phys 2002, 116, 9811
#             [2] Ding et al. Molecular Simulation, 
#                 DOI: 10.1080/08927022.2015.1047367 
#             [3] Weast, Handbook of chemistry and physics, 61st ed, CRC
#                 Boca Raton, FL, 1980

import sys

if len(sys.argv) != 5:
  print "Syntax: python calcHeatCap.py E1 E2 T1 T2"
  sys.exit()

E1 = float(sys.argv[1])
E2 = float(sys.argv[2])
T1 = float(sys.argv[3])
T2 = float(sys.argv[4])

Q = -2.22 # [kcal/mol] quantum contribution @298 K (see ref [1] above)

print "%f kcal/mol" % ( (E2-E1)/(T2-T1) + Q )
