
reference = SLP_R(0.001, 10, 200, 4, [pi, pi], [1, 2], 1);
num_points = 90; % Number of iterations

errors = zeros(num_points/3, 1); % Array to store errors

start = 11;

x_values = start:3:start + 3*(num_points/3 - 1); % Define x-values for the data points

for i = start:3:num_points 
    mat = SLP_R(0.001, 10, i, 4, [pi, pi], [1, 2], 1);
    diff = reference - mat;
    err = norm(diff, 2);
    errors((i-start)/3+1) = err;
end 

% Plotting the error at each iteration
figure;
semilogy(x_values, errors, '-o');
xlabel('Iteration');
ylabel('Error');
title('Error at each Iteration');


function SLP_matrix = SLP_R(w, N_lattice, num_points, N_SLP, alpha, beta, R)
    SLP_size = 2 * N_SLP + 1;
    %num_points = num_points*5;
    SLP_matrix = zeros(SLP_size , SLP_size);            % Initialize matrix
    theta_values = linspace(0, 2*pi, num_points + 1);   % Discretized interval from 0 to 2*pi
    phi_values = linspace(0, 2*pi, num_points + 1);     % Discretized interval from 0 to 2*pi
    theta_values = theta_values(1:end-1);               % prevent counting it twice since mod(2pi)
    phi_values = phi_values(1:end-1);                   % prevent counting it twice since mod(2pi)
    G_matrix = zeros(num_points, num_points);           % Initialize the matrix to store the results

    % Compute the matrix G_{t', t} which is of the size as the
    for i = 1:num_points                        % discretization of the boundary integral
        for j = 1:num_points
            % x = R*(exp(1i*theta_values(i))-exp(1i*phi_values(j))); 
            x_1 = R * (cos(theta_values(i)) - cos(phi_values(j)));
            x_2 = R * (sin(theta_values(i)) - sin(phi_values(j)));
            x = [x_1; x_2];                     % lattice sum function expects a vector
            scale = (2*pi*R)^2/num_points^2;    % Prefactor stems from the quadrature rule.
            G_matrix(i, j) = scale *  lattice_sum(x, w, N_lattice, alpha, beta, R); % the matrix entries are complex
       end
    end
    
    % Compute entries of the SLP matrix
    m_range = - N_SLP:N_SLP;        % Finite Fourier Basis
    n_range = - N_SLP:N_SLP;
    
    for i = 1:SLP_size              % SLP_size = 2 * N_SLP + 1
        m = m_range(i);             % i.e. m = {-N_SLP,...,0,N_SLP}
        for j = 1:SLP_size
            n = n_range(j);

            SLP_matrix(i, j) = exp(-1i * m * theta_values) * (G_matrix * transpose(exp(1i * n * phi_values)));
       end
    end

end


function SLP = makeR(k0, R, alpha, beta, num_points, N_SLP, N_lattice)
    rho_0 = 1000;
    kappa_0 = 1000;
    v0 = sqrt(kappa_0/rho_0);
    w = k0 / v0;                % rescale from k and w
    SLP = SLP_R(w, N_lattice, num_points, N_SLP, alpha, beta, R);
end
 
% Generates the complex greens function f(q) returns a complex scalar
function result = Green(q, x, w, alpha, beta)
    q_alpha = alpha + q;
    
    denominator1 = norm(q_alpha,2)^2 + 2j * dot(beta, q_alpha) - w^2 - norm(beta,2)^2;
    denominator2 = norm(q_alpha,2)^2 - w^2;
    
    numerator = norm(beta,2)^2 - 2j * dot(beta, q_alpha);
    
    result = exp(1j * dot(q_alpha, x)) * (numerator / (denominator1 * denominator2));
end

% Computes the lattice sum G^{\alpha, \beta}(x) returns a scalar
function lattice_sum_value = lattice_sum(x, w, n_lat, alpha, beta, R)
    a = [2*pi, 0]; % Initialize the dual lattice
    b = [0, 2*pi];
    q_0 = [0, 0];  % Center of the lattice
    lattice_sum_value = 0;
    for i = -n_lat:n_lat
        for j = -n_lat:n_lat
            q = q_0 + i * a + j * b;
            lattice_sum_value = lattice_sum_value + Green(q, x, w, alpha, beta);
            % if q ~= 0.1 (test used for debugging)
            %    lattice_sum_value = lattice_sum_value + f(q, x, w, alpha, beta);
            % end
        end
    end    
    lattice_sum_value = 1/(1) * lattice_sum_value; %rescale to volume of the uni cell, (|Y| = 1)
end
