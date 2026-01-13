# Spin_molecular_dynamics_by_magnetic_microwaves

## Sintering Experiments and Theory

Sintering by experiments, Prof. R. Roy et. al., Nature, 399, 688 (1999).  

Theory identifying unpaired 3d electron spins to increase above Curie temperture,  
Selective heating mechanism of magnetic metal oxides by a microwave 
magnetic field, Phys. Rev. B, 79, 104420 (2009).


## Simulation Procedure

Moving spins of Fe(3+), Fe(2+) and O(2-) in cubic cells

Microwaves of giga-Hertz frequency

Dissipated spin molecular dynamics simulation 

  > execution of 1,000,000 steps

  > Metropolis criterion for the next step: accept/reject




## Numerical Code, Parameters and Files

@spin_nucCLD7M6.f03: numerical code

param_spinRL6.h: parameters 

 > number of nodes, total of irons and oxygens, p3m resolution
  
SAI106_config.START1: configuration

 > physical run time, lattice size, number of cells, exchange integrals for Fe and O,
  period of microwave magnetic field, temperature, Curie temprature, etc.
  About 1,000,000 steps are required !

magnetite8.xyz: magnetite in a cubic lattice



