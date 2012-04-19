#!/usr/bin/env python

# Script: calcNeighFrac.py
# Author: Mario Orsi (orsimario at gmail.com)
# Purpose: Reads a LAMMPS ".trj" trajectory (dump) file and, for a
#          user-defined type "i", calculates the fraction of nearest
#          neighbors of the same type "i" in the xy-plane
# Notes: - this script works for lipid bilayers arranged parallel to
#          the xy-plane, and it might need modifications to handle
#          different systems
#        - all and only the nearest 4 neighbors are considered (but
#          it should be easy to modify/generalize the script)
#        - the calculation is done separately for the two monolayers,
#          furthermore it is assumed that the input trajectory files
#          contains z coordinates that centered at 0 (so that the two
#          monolayers can be identify from their positive vs negative
#          coordinates)
# Syntax: calcNeighFrac.py trjFile [upper|lower] type typeNeigh
# Example: calcNeighFrac.py dump1000.trj 2 7
# Reference: de Vries et al, J Phys Chem B 2004, 108, 2454

import sys, string, linecache
from math import sqrt

# Wrap function: wrap a (distance) vector around periodic boundaries
def Wrap( distance, edge ):
  if distance > 0.5*edge:
    distance -= edge
  elif distance < -0.5*edge:
    distance += edge
  return distance

if len(sys.argv) != 5:
  print "Syntax: calcNeighFrac.py trjFile [upper|lower] type typeNeigh"
  sys.exit()

inFileName = sys.argv[1]
monolayer = sys.argv[2]
type = int(sys.argv[3])
typeNeigh = int(sys.argv[4])
inFile = open(inFileName, "r")

if monolayer != "upper" and monolayer != "lower":
  print "Syntax error: either 'upper' or 'lower' layer must be specified"
  sys.exit()

print "Processing file %s ..." % inFileName
line = linecache.getline(inFileName, 4)
words = string.split(line)
nAtoms = int(words[0])
print "Number of \'atoms\': %d" % nAtoms

# get x and y box dimensions
line = linecache.getline(inFileName, 6)
words = string.split(line)
xBox = abs(float(words[0])-float(words[1]))
print "x dimension: %f A" % xBox

line = linecache.getline(inFileName, 7)
words = string.split(line)
yBox = abs(float(words[0])-float(words[1]))
print "y dimension: %f A" % yBox

lines = inFile.readlines()
nAtoms=0; # atom counter
x = [] # x coord 
y = [] # y coord 
z = [] # z coord 
t = [] # type
t1n = [] # type of 1st neighbor
t2n = [] # type of 2nd neighbor
t3n = [] # type of 4th neighbor
t4n = [] # type of 5th neighbor

print "Scanning atoms..."
for line in lines:
  words = string.split(line)
  if len(words) == 9: # (relies on assumption on input file format)
    tAtom=int(words[1]) # type identifier
    xCoord=float(words[3])
    yCoord=float(words[4])
    zCoord=float(words[5])
    # store properties of atoms of the desired type(s)
    if ((tAtom == type) or (tAtom == typeNeigh)):
      nAtoms += 1
      x.append(xCoord)
      y.append(yCoord)
      z.append(zCoord)
      t.append(tAtom)

if monolayer == "upper":
  switch = 1.
elif monolayer == "lower":
  switch = -1.

for i in range( 0, nAtoms ): # loop over stored atoms
  d1i = 999999999.9 # initialize 1st neigh distance
  d2i = 999999999.9 # initialize 2nd neigh distance
  d3i = 999999999.9 # initialize 3rd neigh distance
  d4i = 999999999.9 # initialize 4th neigh distance
  t1n.append(999999999) # initialize 1st neigh type
  t2n.append(999999999) # initialize 2nd neigh type
  t3n.append(999999999) # initialize 3rd neigh type
  t4n.append(999999999) # initialize 4th neigh type
  if z[i]*switch > 0.: # consider only upper or lower layer
    for j in range( 0, nAtoms ): # loop over pairs
      if i!=j: # exclude same-atom pair
        xij = x[i]-x[j] # compute pair x-distance
        yij = y[i]-y[j] # compute pair y-distance
        # perform periodic wraparound if needed:
        xij = Wrap( xij, xBox )
        yij = Wrap( yij, yBox )     
        # compute xy distance:
        dij = sqrt( xij**2 + yij**2 )
        
        # test distance and (re)rank neighbors:
        if dij < d1i: # if j is a new nearest 1st neighbor to i
          d4i = d3i # new 4th neigh dist is old 3rd neigh dist
          t4n[i] = t3n[i] # new 4th neigh is old 3rd neigh
          d3i = d2i # new 3rd neigh dist is old 2nd neigh dist
          t3n[i] = t2n[i] # new 3rd neigh is old 2nd neigh
          d2i = d1i # new 2nd neigh dist is old 1st neigh dist
          t2n[i] = t1n[i] # new 2nd neigh is old 1st neigh
          d1i = dij # new 1st neigh dist
          t1n[i] = t[j] # new 1st neighbor
        elif dij < d2i: # if j is a new nearest 2nd neighbor to i
          d4i = d3i # new 4th neigh dist is old 3rd neigh dist
          t4n[i] = t3n[i] # new 4th neigh is old 3rd neigh
          d3i = d2i # new 3rd neigh dist is old 2nd neigh dist
          t3n[i] = t2n[i] # new 3rd neigh is old 2nd neigh
          d2i = dij # new 2nd neigh dist
          t2n[i] = t[j] # new 2nd neigh
        elif dij < d3i: # if j is a new nearest 3rd neighbor to i
          d4i = d3i # new 4th neigh dist is old 3rd neigh dist
          t4n[i] = t3n[i] # new 4th neigh dist is old 3rd neigh dist
          d3i = dij # new 3rd neigh dist
          t3n[i] = t[j] # new 3rd neigh
        elif dij < d4i: # if j is a new nearest 4th neighbor to i
          d4i = dij # new 4th neigh dist 
          t4n[i] = t[j] # new 4th neigh

nSameNeighs = 0. # init counter for same-type (as "i") neighs 
nDiffNeighs = 0. # init counter for different-type (than "i") neighs 
for i in range( 0, nAtoms ): # loop over stored atoms and sort neighs
  if z[i]*switch > 0: # consider only upper or lower layer
    # 1st (nearest) neigh:
    if t1n[i] == type:
      nSameNeighs += 1
    else:
      nDiffNeighs += 1 
    # 2nd neigh:
    if t2n[i] == type:
      nSameNeighs += 1
    else:
      nDiffNeighs += 1 
    # 3rd neigh:  
    if t3n[i] == type:
      nSameNeighs += 1
    else:
      nDiffNeighs += 1 
    # 4th neigh:
    if t4n[i] == type:
      nSameNeighs += 1
    else:
      nDiffNeighs += 1

fraction = nSameNeighs / ( nSameNeighs + nDiffNeighs )
print "Fraction of %d-%d neighbors = %f" % (type, type, fraction)
