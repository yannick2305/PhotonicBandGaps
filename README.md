# Complex Band Structure

In this computational notebook, we consider the quasiperiodic Helmholtz scattering problem:

$$\hspace{10mm}\Delta v + k^2u = 0, \quad \text{in } Y \setminus \overline{D}.$$  
$$\hspace{3mm}\Delta v + k_i^2u = 0, \quad \text{in } D_i.$$  
$$u\rvert_+ - u\rvert_- =0, \quad \text{on } \partial D.$$  
$$\frac{\partial u}{\partial \nu} \vert_{-} - \delta \frac{\partial u}{\partial \nu} \vert_{+} = 0, \quad \text{on } \partial D.\hspace{7mm}$$  
$$u(x + \ell) = e^{i (\alpha+i \beta) \cdot \ell}u(x), \quad \text{for all } \ell \in \Lambda.\hspace{3mm}$$  



## One Dimensional Resonator Chains:

<p align="center"> <img src="Figures/1dchain.png" alt="1dchain" width="600"/> </p>

## Bandfunctions:

### Subwavelength Regime
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

- `OneD_InterfaceModes.m`

### Non-Hermitian Skinn effect
- `OneD_SkinEffect.m`
- `OneD_RandomGauge.m`
- `OneD_AlgebraicSkin.m`

<p align="center"> <img src="Figures/OneSkin.png" alt="OneSkin" width="300"/> </p>

We consider the non-Hermitian Skinn effect in a variety of settings:
The **Run Files** for one-dimensional chains are named:
- `OneD_(...).m`

## Two Dimensional Resonator Chains:

### Setup:
We consider a two-dimensional infinite screen of circular resonators.

<p align="center"> <img src="Figures/Resonator_screen.png" alt="1dscreen" width="600"/> </p>

### Multipole Expansion:
In the case of spherical resonators, the solutions to the scattering problem are expressed in spherical harmonics. By choosing a finite-order multipole expansion, the single-layer potential admits a matrix representation [Section A.2, 1].

### Quasiperiodic Capacitance Matrix:
The Capacitance Matrix is a computationally efficient way to reduce the scattering problem to a finite eigenvalue problem. The eigenvalues of the quasiperiodic capacitance matrix then parameterise the band functions of the spectral problem [Theorem 3.8, 1].

## Singularities in the Band Structure
- `TwoD_CapacitanceSurface.m`
<p align="center"> <img src="Figures/Bandsurface.png" alt="Bandsurface" width="300"/> </p>

-`TwoD_FieldSol.m`
<p align="center"> <img src="Figures/FIeldsol.png" alt="Fieldsol" width="300"/> </p>

`TwoD_SLP_KernelSurface.m`
<p align="center"> <img src="Figures/SLP_Surface.png" alt="SLP_Surface" width="300"/> </p>

## Defect modes
- `TwoD_ComplexBands.m`
<p align="center"> <img src="Figures/BandDecay.png" alt="BandDecay" width="500"/> </p>


## Runtime and Convergence
- `TwoD_RuntimeSLP.m`
  <p align="center"> <img src="Figures/Runtime.png" alt="Runtime" width="300"/> </p>

- `TwoD_ConvergenceSLP.m`
  <p align="center"> <img src="Figures/Convergence.png" alt="Convergence" width="300"/> </p>
  

## References:
When using the code in the repository, please cite the following two references:


[1] De Bruijn Y and Hiltunen E O, *Complex Band Structure for Subwavelength Evanescent Waves*. Studies in Applied Mathematics (https://doi.org/10.48550/arXiv.2409.14537).

[2] De Bruijn Y and Hiltunen E O, *Complex Band Structure for non-Hermitian spectral problems*.
