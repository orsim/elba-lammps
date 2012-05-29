units		real
atom_style	hybrid angle dipole sphere 
read_data 	data.128lip4232wat
include 	forcefield.elba
velocity	all create 0.0 87287 

neighbor	1 bin
neigh_modify	delay 0

fix		integrate all nve/sphere update dipole
fix 		thermo all langevin 303.15 303.15 1000 9 omega yes zero yes
fix             removeMomentum all momentum 100 linear 1 1 1

thermo		10
dump		trj all custom 50 dump.*.trj id type mol x y z mux muy muz
dump_modify	trj sort id pad 6

timestep	2
run		100

timestep	10
run		100

fix		baro all press/berendsen aniso 1 1 500 modulus 21740

run 		100