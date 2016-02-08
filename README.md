The ELBA-LAMMPS toolkit is a collection of tools aimed at facilitating the
simulation of the ELBA model in LAMMPS. More concretely, ELBA-LAMMPS can be 
described as a LAMMPS add-on package providing users with various 
functionalities to setup, run and analyze molecular dynamics simulations 
of systems modeled with the ELBA force field.

The ELBA model is described in the following papers:

* _Stress testing the ELBA water model_,
  W. Ding, M. Palaiokostas, and M. Orsi,
  [Mol. Simul. 42, 1 (2016)](http://dx.doi.org/10.1080/08927022.2015.1047367).
* _Comparative assessment of the ELBA coarse-grained model for water_,
  M. Orsi,
  [Mol. Phys. 112, 1 (2014)](http://dx.doi.org/10.1080/00268976.2013.844373).
* _Direct Mixing of Atomistic Solutes and Coarse-Grained Water_,
  M. Orsi, W. Ding, and M. Palaiokostas,
  [J. Chem. Theory Comput. 10, 4684 (2014)](http://dx.doi.org/10.1021/ct500065k).
* _Physical properties of mixed bilayers containing lamellar and nonlamellar lipids:
  insights from coarse-grain molecular dynamics simulations_,
  M. Orsi and J. W. Essex,
  [Faraday Discuss. 161, 249 (2013)](http://dx.doi.org/10.1039/c2fd20110k).
* _The ELBA force field for coarse-grain modeling of lipid membranes_,
  M. Orsi and J. W. Essex,
  [PLoS One 6, e28637 (2011)](http://dx.doi.org/10.1371/journal.pone.0028637).

The various folders in the ELBA-LAMMPS toolkit contain usage information 
in the form of README files. 

- examples: various test runs, which can be used to check the LAMMPS build,
	to understand how to run and analyze ELBA systems, and as starting 
	points for new simulations 

- tools: scripts performing various operations, such as analysis of LAMMPS 
	output files to extract properties of interest and file manipulation
	to enable trajectory visualization; check comments at the top of 
	each file to find usage information 

- viz: status files for visualization of ELBA systems with VMD

Mario Orsi (m.orsi@qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
