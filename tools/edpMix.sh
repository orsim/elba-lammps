#!/bin/bash 
numDens2elecDens.py numDensWat.zProfile 10 > edpWat.dat
numDens2elecDens.py numDensChol.zProfile 50 > edpChol.dat
numDens2elecDens.py numDensPhos.zProfile 47 > edpPhos.dat
numDens2elecDens.py numDensGly.zProfile 39 > edpGly.dat
numDens2elecDens.py numDensEst.zProfile 30 > edpEst.dat
numDens2elecDens.py numDensTail.zProfile 23.8 > edpTail.dat
numDens2elecDens.py numDensAmi.zProfile 26 > edpAmi.dat
sumProfile.py edp* > edpSum.dat
