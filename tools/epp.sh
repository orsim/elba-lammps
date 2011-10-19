#!/bin/bash 
mu2Psi.py muzWat.zProfile 4600 > eppWat.dat
mu2Psi.py muzGly.zProfile 4600 > eppGly.dat
mu2Psi.py muzEst.zProfile 4600 > eppEst.dat
qDens2Psi.py qDensHead.zProfile 4600 > eppHead.dat
sumProfile.py epp* > eppSum.dat
