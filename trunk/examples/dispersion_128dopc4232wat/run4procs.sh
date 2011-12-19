#!/bin/bash          
mpirun -np 4 lammps -echo both < in.self
