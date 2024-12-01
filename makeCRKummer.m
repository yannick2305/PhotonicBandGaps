%{
    ------------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Generate complex capacitance matrix using Kummer's method]
    ------------------------------------------------------------------------------
%}

function matC = makeCRKummer(k, R,  beta, N_multi, N_lattice)  
    %
    % Description:
    %   Function to compute the capacitance matrix using Single Layer Potential
    %   and Multipole Expansion and Kummer's transformation.
    
    % --- Generate the complex Single Layer Potential ---
    matSR = makeSKa0(k, R, beta, N_multi, N_lattice);

    % --- Multipole Expansion ---
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
