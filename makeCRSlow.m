%{
    --------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Generate complex Capacitance matrix using for loops]
    --------------------------------------------------------------------
%}

function matC = makeCRSlow(k, R, alpha, beta, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice)  
    %
    % Description:
    %   Function to compute the capacitance matrix using Single Layer Potential
    %   and Multipole Expansion.
    %
    % Input Arguments:
    %   k          - Wave number (scalar).
    %   R          - Resonator radius (scalar).
    %   alpha      - Real- quasimomentum, vector (R^2).
    %   beta       - Complex part of the quasimomentum, vector (R^2).
    %   L1x        - Lattice size in the x-direction (scalar).
    %   L2         - Lattice size in the y-direction (scalar).
    %   d_zeta     - Lattice spacing (scalar).
    %   JHdata     - Data for JH matrix (structure or matrix).
    %   JHijdata   - Data for JHij matrix (structure or matrix).
    %   N=1        - Number resonators inside unit cell (scalar).
    %   N_multi    - Multipole expansion range (scalar).
    %   N_lattice  - Number of lattice points in Green's function (scalar).
    %
    % Output Arguments:
    %   matC       - Capacitance matrix (scalar).
 

    N = 1; % This code only works for one resonator inside the unit cell.

    % -------- Generate the complex Single Layer Potential -----------
    matS  = makeS(k, R, alpha, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    matR  = makeRSlow(k, R, alpha, beta, N_multi, N_lattice);
    matSR = (-matR + matS);  % Minus sign due to different conventions (Habib vs Erik)

    % --- Multipole Expansion -----------------------------------------
    n_range = -N_multi : N_multi;  % Fourier basis representation

    % --- Precompute constants ---
    theta1  = atan2(beta(2), beta(1));
    theta2  = atan2(-beta(2), -beta(1));
    Rb_norm = R * norm(beta);

    % --- Vectorized Fourier coefficient computation ---
    n  = n_range(:);                                    % Ensure n_range is a column vector
    bI = besseli(n, Rb_norm);                           % Compute Bessel functions for all n
    coefficients_P = bI .* exp(-1i * n * theta1);       % Positive Fourier coefficients
    coefficients_M = bI .* exp(-1i * n * theta2);       % Negative Fourier coefficients

    % --- Compute scalar product to define the Capacitance matrix ---
    phi_j = coefficients_P; 
    psi_j = matSR \ phi_j;                              % Solve for psi_j by inverting matSR
    matC  = -2 * pi * R * dot(coefficients_M, psi_j);   % Final capacitance matrix computation

end
