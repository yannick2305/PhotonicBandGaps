%{
    ------------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Generate Complex capacitance without accelerated lattice sums]
    ------------------------------------------------------------------------------
%}



function matC = makeCRIntuitive(k, R, alpha, beta, N_multi, N_lattice)  
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
    %   N_multi    - Multipole expansion range (scalar).
    %   N_lattice  - Number of lattice points in Green's function (scalar).
    %
    % Output Arguments:
    %   matC       - Capacitance matrix (scalar).
 

    % -------- Generate the complex Single Layer Potential -----------
    matSR = -makeRIntuitive(k, R, alpha, beta, N_multi, N_lattice*5);

    % --- Multipole Expansion -----------------------------------------
    n_range = -N_multi : N_multi;  % Fourier basis representation

    % --- Precompute constants ---
    theta1 = atan2(beta(2), beta(1));
    theta2 = atan2(-beta(2), -beta(1));
    Rb_norm = R * norm(beta);

    % --- Vectorized Fourier coefficient computation ---
    n = n_range(:);                                     % Ensure n_range is a column vector
    bI = besseli(n, Rb_norm);                           % Compute Bessel functions for all n
    coefficients_P = bI .* exp(-1i * n * theta1);       % Positive Fourier coefficients
    coefficients_M = bI .* exp(-1i * n * theta2);       % Negative Fourier coefficients

    % --- Compute scalar product to define the Capacitance matrix ---
    phi_j = coefficients_P; 
    psi_j = matSR \ phi_j;                              % Solve for psi_j by inverting matSR
    matC = -2 * pi * R * dot(coefficients_M, psi_j);    % Final capacitance matrix computation

end
