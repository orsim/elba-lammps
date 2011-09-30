units		real
atom_style	hybrid angle dipole sphere 

read_data 	data.in
include 	forcefield.dopc

special_bonds	lj/coul 0.0 1.0 1.0 

velocity	all create 0.0 87287 

neighbor	1 bin
neigh_modify	delay 0

group		lipids type 2 3 4 5 6
group		water type 1

fix		integrate all nve/sphere update dipole

fix 		tempLipids lipids langevin 298.15 298.15 100.0 48279 omega yes
fix 		tempWater water langevin 298.15 298.15 100.0 48279 omega yes

timestep	1
run		50

