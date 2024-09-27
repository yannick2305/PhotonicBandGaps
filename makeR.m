

%% Test the function parameters: (debugging purposes)

 %w = 0.58;             % Set a resonant frequency
 %R = 1;                % Set the radius of the resonator
 %alpha = [pi, pi];     % Fix alpha to the edge of the Brillouin zone
 %beta = [pi, pi];      % Beta will be varied later on in the plot
 %N_SLP = 2;            % The size of the truncated Fourier Basis
 %N_lattice = 20;       % Size of the truncated lattice

 %%Call the makeR function with the provided variables
 %SLP_result = makeR_l(w, R, alpha, beta, N_SLP, N_lattice);

 %%Display or use the result as needed
 %disp(SLP_result);

%% Define the functions:

function SLP = makeR(k0, R, alpha, beta, N_SLP, N_lattice)
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

