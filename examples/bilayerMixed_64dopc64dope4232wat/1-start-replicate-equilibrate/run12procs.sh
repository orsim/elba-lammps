#!/bin/bash          
mpirun -np 12 lammps -echo both < in.start-replicate-equilibrate
