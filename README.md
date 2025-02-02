# Complex Band Structure

## One Dimensional Resonator Chains:

The **Run Files** for one-dimensional chains are named:
- `OneD_(...).m`

### Transfer/Propagator Matrix:
In one-dimensional systems, the Helmholtz scattering problem reduces to an ODE, meaning that its solution can be propagated from initial values. The propagation of an eigenmode over one unit cell is modeled via the transfer matrix.

## Two Dimensional Resonator Chains:
The **Run Files** for two-dimensional chains are named:
- `TwoD_(...).m`

### Setup:
We consider a two-dimensional infinite screen of circular resonators.

### Multipole Expansion:
In the case of spherical resonators, the solutions to the scattering problem are expressed in spherical harmonics. By choosing a finite-order multipole expansion, the single-layer potential admits a matrix representation [Section A.2, 1].

### Quasiperiodic Capacitance Matrix:
The Capacitance Matrix is a computationally efficient way to reduce the scattering problem to a finite eigenvalue problem. The eigenvalues of the quasiperiodic capacitance matrix then parameterise the band functions of the spectral problem [Theorem 3.8, 1].

## References:
When using the code in the repository, please cite the following two references:


[1] De Bruijn Y and Hiltunen E O, *Complex Band Structure for Subwavelength Evanescent Waves*. Studies in Applied Mathematics (https://doi.org/10.48550/arXiv.2409.14537).

[2] De Bruijn Y and Hiltunen E O, *Complex Band Structure for non-Hermitian spectral problems*.
