% compute the single layer potential for a dimer


function SLP = makeR_2(k0, z, R, alpha, beta, N_SLP, N_lattice)
    SLP = SLP_R(k0, z, N_lattice, N_SLP, alpha, beta, R);
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
function [ls1, ls2, ls3] = lattice_sum(m, n, z, w, n_lat, alpha, beta, R)
    a = [2*pi, 0]; % Initialize the dual lattice
    b = [0, 2*pi];
    sum_1 = 0;
    sum_2 = 0;
    sum_3 = 0;
    for i = -n_lat:n_lat
        for j = -n_lat:n_lat
            q = i * a + j * b;
            alpha_q = alpha + q;
            psi = atan2(alpha_q(2), alpha_q(1));
            scale = 2*pi*R * exp(1i * psi * (n - m)) * 1i^(m - n) * (-1)^n;
            arg = R*sqrt(alpha_q(1)^2 + alpha_q(2)^2);
            res = scale * besselj(m, arg) * besselj(n, arg)* Green(q, w, alpha, beta);
      
            sum_1 = sum_1 + res;
            sum_2 = sum_2 + res*exp(-1i*dot(alpha_q,z));
            sum_3 = sum_3 + res*exp(1i*dot(alpha_q,z));
        end
    end    
    ls1 = 1/(1) * sum_1; %rescale to volume of the unit cell, (|Y| = 1)
    ls2 = 1/(1) * sum_2; %rescale to volume of the unit cell, (|Y| = 1)
    ls3 = 1/(1) * sum_3; %rescale to volume of the unit cell, (|Y| = 1)
end

% Perform the boundary integral for SLP, returns of form Mat(2N+1)
function SLP_matrix = SLP_R(w, z, N_lattice, N_SLP, alpha, beta, R)
    SLP_size =  (2 * N_SLP + 1);
 
    SLP_matrix_1 = zeros(SLP_size , SLP_size);            % Initialize matrix
    SLP_matrix_2 = zeros(SLP_size , SLP_size);
    SLP_matrix_3 = zeros(SLP_size , SLP_size);


    % Compute entries of the SLP matrix
    m_range = - N_SLP:N_SLP;        % Finite Fourier Basis
    n_range = - N_SLP:N_SLP;
    
    for i = 1:SLP_size              % SLP_size = 2 * N_SLP + 1
        m = m_range(i);             % i.e. m = {-N_SLP,...,0,...,N_SLP}
        for j = 1:SLP_size
            n = n_range(j);
            [l1,l2,l3] = lattice_sum(m, n, z, w, N_lattice, alpha, beta, R);
            SLP_matrix_1(i, j) = l1;
            SLP_matrix_2(i, j) = l2;
            SLP_matrix_3(i, j) = l3;
       end
    end

    %SLP_matrix_3 = transpose(SLP_matrix_2);
    

    SLP_matrix = [SLP_matrix_1, SLP_matrix_2 ; SLP_matrix_3, SLP_matrix_1];



end

