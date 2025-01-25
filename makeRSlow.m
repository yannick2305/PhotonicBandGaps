%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [January 2025]
    Description:  [Non-vectorised computation of the SLP]
    ----------------------------------------------------------------------
%}

function SLP = makeRSlow(k0, R, alpha, beta, N_SLP, N_lattice)
    rho_0 = 1000;
    kappa_0 = 1000;
    v0 = sqrt(kappa_0/rho_0);
    w = k0 / v0;                % rescale from k and w (Here k = w)
    SLP = SLP_R(w, N_lattice, N_SLP, alpha, beta, R);
    SLP = SLP_R(k0, N_lattice, N_SLP, alpha, beta, R);
end

% Generates part of the complex Green's function
function result = Green(q, w, alpha, beta)
    q_alpha = alpha + q;

    denominator1 = norm(q_alpha,2)^2 + 2j * dot(beta, q_alpha) - w^2 - norm(beta,2)^2;
    denominator2 = norm(q_alpha,2)^2 - w^2;

    numerator = norm(beta,2)^2 - 2j * dot(beta, q_alpha);

    result = (numerator / (denominator1 * denominator2));
end

% Computes the Single Layer potential at a given entry
function lattice_sum_value = lattice_sum(m, n, w, n_lat, alpha, beta, R)
    a = [2*pi, 0]; % Initialize the dual lattice
    b = [0, 2*pi];
    lattice_sum_value = 0;
    for i = -n_lat:n_lat
        for j = -n_lat:n_lat
            q = i * a + j * b;
            alpha_q = alpha + q;
            psi = atan2(alpha_q(1), alpha_q(2));
            psi = atan2(alpha_q(2), alpha_q(1));
            scale = 2*pi*R * exp(1i * psi * (n - m)) * 1i^(m - n) * (-1)^n;
            arg = R*sqrt(alpha_q(1)^2 + alpha_q(2)^2);
            res = scale * besselj(m, arg) * besselj(n, arg)* Green(q, w, alpha, beta);

            lattice_sum_value = lattice_sum_value + res;

        end
    end    
    lattice_sum_value = 1/(1) * lattice_sum_value; %rescale to volume of the unit cell, (|Y| = 1)
end

% Perform the boundary integral for SLP, returns of form Mat(2N+1)
function SLP_matrix = SLP_R(w, N_lattice, N_SLP, alpha, beta, R)
    SLP_size = 2 * N_SLP + 1;

    SLP_matrix = zeros(SLP_size , SLP_size);            % Initialize matrix

    % Compute entries of the SLP matrix
    m_range = - N_SLP:N_SLP;        % Finite Fourier Basis
    n_range = - N_SLP:N_SLP;

    for i = 1:SLP_size              % SLP_size = 2 * N_SLP + 1
        m = m_range(i);             % i.e. m = {-N_SLP,...,0,...,N_SLP}
        for j = 1:SLP_size
            n = n_range(j);
            SLP_matrix(i, j) = lattice_sum(m, n, w, N_lattice, alpha, beta, R);
       end
    end

end
