# Spin_molecular_dynamics_by_magnetic_microwaves

## Sintering Experiments and Theory

Sintering by initial experiments (R. Roy et. al., Nature, 1999), Ref.1.  

Theory of showing unpaired 3d electron spins of Fe(3+) and Fe(2+) to increase 
above the Curie temperature (M. Tanaka et al., Phys. Rev. B, 2009), Ref.2.  


## Simulation Procedure

Electron spins of Fe(3+), Fe(2+) and O(2-) in cubic cells

Microwaves of giga-Hertz frequency, like 2.5 GHz

Dissipated spin molecular dynamics simulation 

  > Execution of a few 1,000,000 steps of $ Delta t= 0.001 $ ps

  > Metropolis criterion to accept/reject for the next steps


## Numerical Code, Parameters and Files

1) @spin_nucCLD7M6a.f03: numerical code

2) param_spinRL6.h: parameters 

 > number of nodes, total of irons and oxygens, p3m resolution
  
3) SAI106_config.START1: configuration

 > physical run time, lattice sizes, number of cells, exchange integrals for Fe and O,
  period of microwave magnetic field, temperature, Curie temprature, etc.
 About 3,000,000 steps are required !

4) Magnetite in a cubic lattice: magnetite8.xyz for initialization

## References

1. R. Roy, D. Agrawal, J. Cheng, and S. Gedevanishvili, Full sintering of powdered-metal bodies
in a magnetic field, Nature, 399, 668 (1999).

3. M. Tanaka, H. Kono, and K. Maruyama, Selective heating mechanism of magnetic
metal oxides by a microwave magnetic field, Phys. Rev. B, 79, 104420 (2009).

