%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [January 2025]
    Description:  [Non-vectorised computation of the SLP]
    ----------------------------------------------------------------------
%}


% --- Perform the boundary integral for SLP, returns of form Mat(2N+1) ---
function SLP_matrix = makeRSlow(w, R, alpha, beta, N_SLP, N_lattice)
    % Initialize matrix
    SLP_size = 2 * N_SLP + 1;
    SLP_matrix = zeros(SLP_size , SLP_size);   
    
    % Compute entries of the SLP matrix
    m_range = - N_SLP:N_SLP;        
    n_range = - N_SLP:N_SLP;
    
    for i = 1:SLP_size      
        m = m_range(i);       
        for j = 1:SLP_size
            n = n_range(j);
            SLP_matrix(i, j) = lattice_sum(m, n, w, N_lattice, alpha, beta, R);
       end
    end

end


% --- Generates part of the complex Green's function ---
function result = Green(q, w, alpha, beta)
    % precompute some variables
    q_alpha          = alpha + q;
    q_alpha_norm     = norm(q_alpha, 2); 
    beta_norm_sq     = norm(beta, 2)^2;  
    q_alpha_dot_beta = dot(beta, q_alpha);  
    
    % Compute denominators and numerator
    denominator1 = q_alpha_norm^2 + 2j * q_alpha_dot_beta - w^2 - beta_norm_sq;
    denominator2 = q_alpha_norm^2 - w^2;
    numerator    = beta_norm_sq - 2j * q_alpha_dot_beta;
    
    result = (numerator / (denominator1 * denominator2));
end


% --- Computes the lattice sum at a given entry ---
function lattice_sum_value = lattice_sum(m, n, w, n_lat, alpha, beta, R)
    % Initialize the dual lattice 
    a = [2*pi, 0]; 
    b = [0, 2*pi];

    lattice_sum_value = 0;
    psi_mn = (-1)^n * 1i^(m + n); % Precompute this term

    for i = -n_lat:n_lat 
        for j = -n_lat:n_lat
            q = i * a + j * b;

            alpha_q = alpha + q;
            psi     = atan2(alpha_q(2), alpha_q(1));
            scale   = 2 * pi * R * exp(1i * psi * (n - m)) * psi_mn;
            arg     = R * sqrt(alpha_q(1)^2 + alpha_q(2)^2);
            res     = scale * besselj(m, arg) * besselj(n, arg) * Green(q, w, alpha, beta);
      
            lattice_sum_value = lattice_sum_value + res;

        end
    end    
    lattice_sum_value = 1/(1) * lattice_sum_value; %rescale to volume of the unit cell, (|Y| = 1)
end
