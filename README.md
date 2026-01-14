# Spin_molecular_dynamics_by_magnetic_microwaves

## Sintering Experiments and Theory

Sintering by experiments (R. Roy et. al., Nature, 1999).  

Theory identifying unpaired 3d electron spins to increase above 
Curie temperture (M. Tanaka et al., J. Chem. Phys, 2009).  


## Simulation Procedure

Electron spins moving of Fe(3+), Fe(2+) and O(2-) in cubic cells

Microwaves of giga-Hertz frequency

Dissipated spin molecular dynamics simulation 

  > execution of about 10,000 steps of $ Delta t= 0.001 $ ps

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

## References

1. R. Roy, D. Agrawal, J. Cheng, and S. Gedevanishvili, Full sintering of powdered-metal bodies
in a magnetic field, Nature, 399, 668 (1999).

3. M. Tanaka, H. Kono, and K. Maruyama, Selective heating mechanism of magnetic
metal oxides by a microwave magnetic field, Phys. Rev. B, 79, 104420 (2009).

