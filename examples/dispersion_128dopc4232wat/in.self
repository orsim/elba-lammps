units		real
atom_style	hybrid angle dipole sphere 
read_data 	data.128lip4232wat
include 	forcefield.elba
velocity	all create 0.0 87287 

neighbor	1 bin
neigh_modify	delay 0

group		lip type 2 3 4 5 6
group		wat type 1

fix		integrate all nve/sphere update dipole
fix 		thermoLip lip langevin 303.15 303.15 100.0 48279 omega yes
fix 		thermoWat wat langevin 303.15 303.15 100.0 48279 omega yes
fix             removeMom all momentum 1 linear 1 1 1

thermo		10
dump		trj all custom 50 dump.trj id type mol x y z mux muy muz
dump_modify	trj sort id

timestep	2
run		100

timestep	10
run		100

fix		baro all press/berendsen aniso 1 1 500 modulus 21740

run 		100