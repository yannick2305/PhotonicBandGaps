% compute the single layer potential for a dimer


function SLP = makeR_2(k0, R, alpha, beta, N_SLP, N_lattice)
    rho_0 = 1000;
    kappa_0 = 1000;
    v0 = sqrt(kappa_0/rho_0);
    w = k0 / v0;                % rescale from k and w (Here k = w)
    SLP = SLP_R(w, N_lattice, N_SLP, alpha, beta, R);
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
function lattice_sum_value_1 = lattice_sum_1(m, n, w, n_lat, alpha, beta, R)
    a = [2*pi, 0]; % Initialize the dual lattice
    b = [0, 2*pi];
    lattice_sum_value = 0;
    for i = -n_lat:n_lat
        for j = -n_lat:n_lat
            q = i * a + j * b;
            alpha_q = alpha + q;
            psi = atan2(alpha_q(1), alpha_q(2));
            scale = 2*pi*R * exp(1i * psi * (n - m)) * 1i^(m - n) * (-1)^n;
            arg = R*sqrt(alpha_q(1)^2 + alpha_q(2)^2);
            res = scale * besselj(m, arg) * besselj(n, arg)* Green(q, w, alpha, beta);
      
            lattice_sum_value = lattice_sum_value + res;

        end
    end    
    lattice_sum_value_1 = 1/(1) * lattice_sum_value; %rescale to volume of the unit cell, (|Y| = 1)
end

function lattice_sum_value_2 = lattice_sum_2(m, n, w, n_lat, alpha, beta, R)
    a = [2*pi, 0]; % Initialize the dual lattice
    b = [0, 2*pi];
    lattice_sum_value = 0;
    % boundary of second resonator
    dis = 0.25; % distance between the resonators;

    for i = -n_lat:n_lat
        for j = -n_lat:n_lat
            q = i * a + j * b;
            alpha_q = alpha + q;
            psi = atan2(alpha_q(1), alpha_q(2));
            % Define the integrand function
            integrand = @(theta, phi) exp(-1i*m*theta) .* exp(1i * exp(1i * psi) .* (R * exp(1i * theta) - (dis * cos(phi) + sqrt(R^2 - dis^2 * sin(phi).^2) .* exp(1i * phi)))) .* exp(1i*n*phi);

            % Compute the double integral using integral2 with options
            % (NEEDS FIXING)
            integ = integral2(integrand, 0, 2*pi, 0, 2*pi, 'RelTol',1e-2,'AbsTol',1e-2);

            scale = R^2;
            res = scale * integ * Green(q, w, alpha, beta);
      
            lattice_sum_value = lattice_sum_value + res;

        end
    end    
    lattice_sum_value_2 = 1/(1) * lattice_sum_value; %rescale to volume of the unit cell, (|Y| = 1)
end



% Perform the boundary integral for SLP, returns of form Mat(2N+1)
function SLP_matrix = SLP_R(w, N_lattice, N_SLP, alpha, beta, R)
    SLP_size =  (2 * N_SLP + 1);
 
    SLP_matrix = zeros(2*(2 * N_SLP + 1) , 2*(2 * N_SLP + 1));
    SLP_matrix_1 = zeros(SLP_size , SLP_size);            % Initialize matrix
    SLP_matrix_2 = zeros(SLP_size , SLP_size);


    % Compute entries of the SLP matrix
    m_range = - N_SLP:N_SLP;        % Finite Fourier Basis
    n_range = - N_SLP:N_SLP;
    
    for i = 1:SLP_size              % SLP_size = 2 * N_SLP + 1
        m = m_range(i);             % i.e. m = {-N_SLP,...,0,...,N_SLP}
        for j = 1:SLP_size
            n = n_range(j);
            SLP_matrix_1(i, j) = lattice_sum_1(m, n, w, N_lattice, alpha, beta, R);
       end
    end

    for i = 1:SLP_size              % SLP_size = 2 * N_SLP + 1
        m = m_range(i);             % i.e. m = {-N_SLP,...,0,...,N_SLP}
        for j = 1:SLP_size
            n = n_range(j);
            SLP_matrix_2(i, j) = lattice_sum_2(m, n, w, N_lattice, alpha, beta, R);
       end
    end
    
    SLP_matrix_3 = transpose(SLP_matrix_2);
    

    SLP_matrix = [SLP_matrix_1, SLP_matrix_2 ; SLP_matrix_3, SLP_matrix_1];



end

