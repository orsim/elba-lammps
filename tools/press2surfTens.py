#!/usr/bin/env python

# Script: press2surfTens.py
# Purpose: Calculates the surface tension using the standard formula:
#          0.5*Lz*(Pzz-0.5(Pxx+Pyy))
# Syntax: press2surfTens.py Lz Pxx Pyy Pzz 
# Notes: Pxx, Pyy and Pzz are provided in units of atm (LAMMPS
#        convention for 'units real' style) and Lz in Angstrom
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)


import sys, string, linecache
from math import sqrt

if len(sys.argv) != 5:
  print "Syntax: press2surfTens.py Lz Pxx Pyy Pzz"
  sys.exit()

atm_in_Pa = 101325 # note: 1 Pa = 1 N/m^2
A_in_m = 1e-10 # Angstrom in meter
N_in_mN = 1e3 # Newton in milliNewton

Lz = float(sys.argv[1]) * A_in_m
Pxx = float(sys.argv[2]) * atm_in_Pa 
Pyy = float(sys.argv[3]) * atm_in_Pa 
Pzz = float(sys.argv[4]) * atm_in_Pa

surfTens = 0.5 * Lz * ( Pzz - 0.5*(Pxx+Pyy) ) # [N/m]
surfTens = surfTens * N_in_mN # [mN/m]

print "Surface tension: % .6f mN/m" % surfTens
