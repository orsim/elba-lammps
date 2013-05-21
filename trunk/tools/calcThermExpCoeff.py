#!/usr/bin/env python

# Script: calcThermExpCoeff.py
# Author: Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Calculates thermal expansion coefficient (alpha) [1]
# Syntax: calcThermExpCoeff.py rho1 rho2 T1 T2
#         rho1 = density at T1
#         rho2 = density at T2
#         T1 = temperature [K]
#         T2 = temperature [K]
# Example: calcThermExpCoeff.py 1.002 0.996 298 318
# Note: Experimental value for water is 2.0 * 10^(-4) / K [2]
# References: [1] Glattli et al, J Chem Phys 2002, 116, 9811
#             [2] Water: A Comprehensive Treatise, F. Frank, Ed;
#                 Plenum, New York, 1972

import sys, math

if len(sys.argv) != 5:
  print "Syntax: calcThermExpCoeff.py rho1 rho2 T1 T2"
  sys.exit()

rho1 = float(sys.argv[1])
rho2 = float(sys.argv[2])
T1 = float(sys.argv[3])
T2 = float(sys.argv[4])

print "%f E-4 / K" % ( - 10**4 * math.log(rho2/rho1) / (T2-T1) )
