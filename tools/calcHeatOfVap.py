#!/usr/bin/env python

# Script: calcHeatOfVap.py
# Author: Mario Orsi
# Purpose: Calculates (molar) heat of vaporization from intermolecular
#          potential energy per mole U and temperature T
# Syntax: python calcHeatOfVap.py U T
#         U = intermolecular potential energy per mole [kcal/mol]
#         T = temperature [K]
# Example: python calcHeatOfVap.py -10 298
# Notes: - Benchmark values for the heat of vaporization of water at 298 K:
#        i) 11.0 kcal/mol from 'classical calculation' [1,2]
#        ii) 10.52 kcal/mol from experiment [3]
#        - Some workers include quantum corrections [4], yet others
#          do not [5]
# References: [1] Guillot, J Mol Liq 101, 219 (2002)
#             [2] Guillot & Guissani, J Chem Phys 114, 6720 (2001)
#             [3] Haar, Gallagher & Kell, NBS/NRC Steam tables,
#                 Hemisphere, Washington DC, 1984
#             [4] Glattli et al, J Chem Phys 2002, 116, 9811
#             [5] Wu et al, J Chem Phys 2006, 124, 024503
#             [6] M. Orsi, Mol. Phys., 112, 1566-1576 (2014)

import sys

if len(sys.argv) != 3:
  print "Syntax: python calcHeatOfVap.py U T"
  sys.exit()

U = float(sys.argv[1])
T = float(sys.argv[2])

R = 1.987e-3 # gas constant in kcal/(mol K)
Q = -0.055 # quantum contribution in kcal/mol (see [4] above)

#print "%f kcal/mol" % ( -U + R*T + Q )
print "%f kcal/mol" % ( -U + R*T )
