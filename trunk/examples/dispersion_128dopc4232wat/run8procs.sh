#!/bin/bash          
mpirun -np 8 lammps -echo both < in.self
