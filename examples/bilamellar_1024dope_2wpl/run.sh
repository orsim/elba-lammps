#!/bin/bash          
mpirun -np 4 lammps -echo both -log log.30sep2011.4procs < in.elba
