#!/usr/bin/env python

# Script: calcNeighFrac.py
# Author: Mario Orsi (orsimario at gmail.com)
# Purpose: Reads a LAMMPS ".trj" trajectory (dump) file and, for a
#          user-defined type "i", calculates the fraction of nearest
#          neighbors of the same type "i" in the xy-plane
# Note: this script works for lipid bilayers arranged parallel to the
#       xy-plane, and it would need modifications to handle different
#       systems
# Syntax: calcNeighFrac.py trjFile type typeNeigh
# Example: calcNeighFrac.py dump.trj 2 4
# Reference: de Vries et al, J Phys Chem B 2004, 108, 2454

import sys, string, linecache
from math import sqrt

if len(sys.argv) != 4:
  print "Syntax: calcNeighFrac.py trjFile type typeNeigh"
  sys.exit()

inFileName = sys.argv[1]
type = int(sys.argv[2])
typeNeigh = int(sys.argv[3])
inFile = open(inFileName, "r")

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

x = [] # x coord of atom
y = [] # y coord of atom
t = [] # atom type
t1n = [] # type of 1st neighbor
t2n = [] # type of 2nd neighbor
t3n = [] # type of 4th neighbor
t4n = [] # type of 5th neighbor

print "Scanning atoms..."
for line in lines:
  words = string.split(line)
  if len(words) == 9: # lipids+water system
    n=int(words[0]) # atom identifier
    tAtom=int(words[1]) # type identifier
    m=int(words[2]) # molecule identifier
    xCoord=float(words[3])
    yCoord=float(words[4])
    if ((tAtom == type) or (tAtom == typeNeigh)):
      nAtoms = nAtoms + 1
      x.append(xCoord)
      y.append(yCoord)
      t.append(tAtom)
for i in range( 0, nAtoms ):
  d1i = d2i = d3i = d4i = 999999999.9 # initialize distances
  t1n.append(999999999) # initialize neigh type
  t2n.append(999999999) # initialize neigh type
  t3n.append(999999999) # initialize neigh type
  t4n.append(999999999) # initialize neigh type
  for j in range( 0, nAtoms ):
    if i!=j:
      xij = x[i]-x[j]
      yij = y[i]-y[j]
      # wrap distances over periodic boundaries:
      if xij > 0.5*xBox:
        xij = xij - xBox
      elif xij < -0.5*xBox:
        xij = xij + xBox
      if yij > 0.5*yBox:
        yij = yij - yBox
      elif yij < -0.5*yBox:
        yij = yij + yBox
      
      # compute xy distance:
      dij = sqrt( xij**2 + yij**2 )

      # test distance
      if dij < d1i: # if j is a new nearest 1st neighbor to i
        d4i = d3i # new 4th neigh dist is old 3rd neigh dist
        t4n[i] = t3n[i] # update type
        d3i = d2i # new 3rd neigh dist is old 2nd neigh dist
        t3n[i] = t2n[i] # update type
        d2i = d1i # new 2nd neigh dist is old 1st neigh dist
        t2n[i] = t1n[i] # update type
        d1i = dij # new 1st neigh dist
        t1n[i] = t[j] # store type
      elif dij < d2i: # if j is a new nearest 2nd neighbor to i
        d4i = d3i # new 4th neigh dist is old 3rd neigh dist
        t4n[i] = t3n[i] # update type
        d3i = d2i # new 3rd neigh dist is old 2nd neigh dist
        t3n[i] = t2n[i] # update type
        d2i = dij # new 2nd neigh dist
        t2n[i] = t[j] # update type
      elif dij < d3i: # if j is a new nearest 3rd neighbor to i
        d4i = d3i # new 4th neigh dist is old 3rd neigh dist
        t4n[i] = t3n[i] # update type
        d3i = dij # new 3rd neigh dist
        t3n[i] = t[j] # update type
      elif dij < d4i: # if j is a new nearest 4th neighbor to i
        d4i = dij # new 4th neigh dist 
        t4n[i] = t[j] # update type

nSameNeighs = 0. # number of same-type neighbors (with respect to "i")
nDiffNeighs = 0. # number of different-type neighbors (with respect to "i")
for i in range( 0, nAtoms ):
  if t1n[i] == type:
    nSameNeighs += 1
  else:
    nDiffNeighs += 1 
  if t2n[i] == type:
    nSameNeighs += 1
  else:
    nDiffNeighs += 1 
  if t3n[i] == type:
    nSameNeighs += 1
  else:
    nDiffNeighs += 1 
  if t4n[i] == type:
    nSameNeighs += 1
  else:
    nDiffNeighs += 1
#print nSameNeighs,nDiffNeighs
fraction = nSameNeighs / (nSameNeighs+nDiffNeighs)
#fraction = nSameNeighs / (4.*nAtoms)
print "Fraction of %d-%d neighbors = %f" % (type, type, fraction)
