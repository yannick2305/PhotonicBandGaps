# Complex Band Structure

In this computational notebook, we consider the quasiperiodic Helmholtz scattering problem:

<p align="center"> <img src="Figures/HelmholtzPDE.png" alt="HelmholtzPDE" width="300"/> </p>



## One Dimensional Resonator Chains:

<p align="center"> <img src="Figures/1dchain.png" alt="1dchain" width="600"/> </p>

## Bandfunctions:

### Subwavelength Regime

#### Quasiperiodic Capacitance
In the subwavelength regime the complex quasiperiodic capacitance matrix allows us to find explicit formulas for the band and gap functions.
- `OneD_Monomer_Band.m`
<p align="center"> <img src="Figures/BandMonomer.png" alt="BandMonomer" width="300"/> </p>

- `OneD_Dimer_Band.m`
 <p align="center"> <img src="Figures/BandDimer.png" alt="BandDimer" width="300"/> </p>

### General Regime

#### Transfer/Propagator Matrix:
In one-dimensional systems, the Helmholtz scattering problem reduces to an ODE, meaning that its solution can be propagated from initial values. The propagation of an eigenmode over one unit cell is modeled via the transfer matrix.

- `OneD_General_Band.m`
  <p align="center"> <img src="Figures/BandGeneral.png" alt="BandGeneral" width="300"/> </p>

## Localisation Effects:

### Exponentially localised Interface modes

The defect supports an eigenfrequency in the band gap of the dimer. The compelex quasimomentum associated to this gap frequency accurately preicts the decay length of the defect eigenmode.

- `OneD_InterfaceModes.m`
 <p align="center"> <img src="Figures/InterfaceBand.png" alt="InterfaceBand" width="500"/> </p>
 <p align="center"> <img src="Figures/InterfaceDecay.png" alt="InterfaceDecay" width="300"/> </p>

### Non-Hermitian Skinn effect
In some non-Hermitian systems such as in the skin effect, the energy leakage can be factored out and the spectral problem can be reformulated as a gap problem, where the complex quasimomentum accurately predicts the decay lenghth of the eigenmodes.

- `OneD_SkinEffect.m`
- `OneD_RandomGauge.m`
- `OneD_AlgebraicSkin.m`

<p align="center"> <img src="Figures/OneSkin.png" alt="OneSkin" width="300"/> </p>

## Two Dimensional Resonator Chains:

### Setup:
We consider a two-dimensional infinite screen of circular resonators.

<p align="center"> <img src="Figures/Resonator_screen.png" alt="1dscreen" width="400"/> </p>

### Multipole Expansion:
In the case of spherical resonators, the solutions to the scattering problem are expressed in spherical harmonics. By choosing a finite-order multipole expansion, the single-layer potential admits a matrix representation [Section A.2, 1].

### Quasiperiodic Capacitance Matrix:
The Capacitance Matrix is a computationally efficient way to reduce the scattering problem to a finite eigenvalue problem. The eigenvalues of the quasiperiodic capacitance matrix then parameterise the band functions of the spectral problem [Theorem 3.8, 1].

## Singularities in the Band Structure

At certain quasiperiodicities the single layer potential fails to be invertible, leading to singularities in the bandfunctions.
- `TwoD_CapacitanceSurface.m`
<p align="center"> <img src="Figures/Bandsurface.png" alt="Bandsurface" width="300"/> </p>

-`TwoD_SLP_KernelSurface.m`
<p align="center"> <img src="Figures/SLP_Surface.png" alt="SLP_Surface" width="300"/> </p>

We compute the field solution poised at a parameter valued point such that the single layer potential fails to e invertible.


-`TwoD_FieldSol.m`
<p align="center"> <img src="Figures/FIeldsol.png" alt="Fieldsol" width="300"/> </p>


## Defect modes

A defected resonator lattice supports eigenfrequencies inside the bandgap. The localisation strength of the defect eigenmode is accurately predicted by the complex band structure.
- `TwoD_ComplexBands.m`
<p align="center"> <img src="Figures/BandDecay.png" alt="BandDecay" width="500"/> </p>


## Runtime and Convergence

We consuct a runtime and convergence analysis, with respect to the truncation size of the lattice sum defining the single layer potential.
- `TwoD_RuntimeSLP.m`
  <p align="center"> <img src="Figures/Runtime.png" alt="Runtime" width="300"/> </p>

- `TwoD_ConvergenceSLP.m`
  <p align="center"> <img src="Figures/Convergence.png" alt="Convergence" width="300"/> </p>
  

## References:
When using the code in the repository, please cite the following two references:


[1] De Bruijn Y and Hiltunen E O, *Complex Band Structure for Subwavelength Evanescent Waves*. Studies in Applied Mathematics (https://doi.org/10.48550/arXiv.2409.14537).

[2] De Bruijn Y and Hiltunen E O, *Complex Band Structure for non-Hermitian spectral problems*.
