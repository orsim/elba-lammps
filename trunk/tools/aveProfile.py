#!/usr/bin/env python

# Script: aveProfile.py
# Author: Julien Michel (www.julienmichel.net)
# Purpose: Reads two or more profiles and calculates average
# Syntax: aveProfile.py fileName1 fileName2 ... fileNameN
# Example: aveProfile.py edp1ns.dat edp2ns.dat edp3ns.dat > edpAve.dat

import sys

averages = []

resized = False
for file in sys.argv[1:]:
    stream = open(file,'r')
    buffer = stream.readlines()
    stream.close()
    if not resized:
        #check = buffer[0].split()
        check = len(buffer)
        for x in range(0,check):
            averages.append([])
        resized = True
    #print averages
    for y in range(0,len(buffer)):
        elems = buffer[y].split()
        nprops = len(elems)
        for x in range(0,len(elems)):
            #print elems[x]
            averages[y].append(float(elems[x]))
# Now collapse
for line in averages:
    avgs = []
    for x in range(0,nprops):
        prop = 0
        y = x
        count = 0
        while (y < len(line)):
            prop += line[y]
            count +=1
            #print y,prop
            y += nprops
        prop /= count
        #print prop
        avgs.append(prop)
    #print line
    str = " "
    for val in avgs:
        str += "%f " % val
    print str
    #sys.exit(-1)
