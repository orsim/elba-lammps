#!/bin/bash 
numDens2dat.py numDensWat.zProfile > ndpWat.dat
numDens2dat.py numDensChol.zProfile > ndpChol.dat
numDens2dat.py numDensPhos.zProfile > ndpPhos.dat
numDens2dat.py numDensGly.zProfile > ndpGly.dat
numDens2dat.py numDensEst.zProfile > ndpEst.dat
numDens2dat.py numDensTail.zProfile > ndpTail.dat
numDens2dat.py numDensAmi.zProfile > ndpAmi.dat
sumProfile.py ndp* > ndpSum.dat
