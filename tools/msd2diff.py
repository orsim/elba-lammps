#!/usr/bin/env python

# Script: msd2diff.py
# Author: Mario Orsi (orsimario at gmail.com, www.soton.ac.uk/~orsi)
# Purpose: Reads a file containing timestep & mean-squared displacement,
#          and converts to corresponding diffusion coefficient
# Syntax: msd2diff.py inputFile timeStep nDimensions
#         timeStep = MD integration time step in fs
#         nDimensions = # of dimensions to consider for diffusion coeff,
#                       typically: 3 for 3-dimensional diff (e.g., bulk
#                       water), 2 for 2-dimensional diff (e.g., lipid
#                       'lateral' motion inside a bilayer)
# Notes: - input msd data are in Angstrom^2
# Examples: - msd2diff.py wat.msd 3 > wat.diff
#           - msd2diff.py lip.msd 2 > lip.diff
# References: - Orsi & Essex, PLoS ONE, submitted
#             - Orsi et al, J Phys Condens Matter 22, 155106 (2010),
#               section 5.5
#             - Rapaport, The Art of Molecular Dynamics Simulation
#               (2004), 2nd ed, p.122
#             - Xiang, J Phys Chem B 103, 385 (1999)

import sys,os,string

if len(sys.argv) != 4:
  print "Syntax: msd2diff.py inputFile timeStep nDimensions"
  sys.exit()

inFileName = sys.argv[1]
dt = float(sys.argv[2])
nDims = int(sys.argv[3])

if (nDims!=1 and nDims!=2 and nDims!=3):
  print "Error: nDimensions must be either 1 or 2 or 3"
  print "Syntax: msd2diff.py inputFile timeStep nDimensions"
  sys.exit()

inFile = open(inFileName, "r")
lines = inFile.readlines()
inFile.close()

# read input data:
time = []
msd = []
for line in lines:
    if line[0] != '#': # ignore comments
        words = string.split(line)
        time.append( dt * float(words[0]) )
        msd.append(float(words[1]))

# unit conversion factors:
fs_in_s = 1e-15 # 1 femtosecond = 10^(-15) second
fs_in_mus = 1e-9 # 1 femtosecond = 10^(-6) microsecond
fs_in_ns = 1e-6 # 1 femtosecond = 10^(-6) nanosecond
A2_in_m2 = 1e-20 # 1 Angstrom^2 = 10^(-20) m^2
A2_in_cm2 = 1e-16 # 1 Angstrom^2 = 10^(-16) cm^2
A2_in_nm2 = 1e-2 # 1 Angstrom^2 = 10^(-2) nm^2

# output diff coeff:
print '# Diffusion coefficient D as a function of measurement time'
print '# Measurement time [ns] | D [cm^2/s] | D [nm^2/mus]'
print '0.0 0.0 0.0' # set D=0 for t_measurement=0
for i in range( 1, len(msd) ):
  print ('%g %g %g %g' % (time[i]*fs_in_ns,
                       msd[i]*A2_in_m2/(2*nDims*time[i]*fs_in_s),
                       msd[i]*A2_in_cm2/(2*nDims*time[i]*fs_in_s),
                       msd[i]*A2_in_nm2/(2*nDims*time[i]*fs_in_mus)))
