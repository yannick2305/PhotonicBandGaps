%{
    ------------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Generate complex single layer potential using Kummer's method]
    ------------------------------------------------------------------------------
%}

function SLP_matrix = makeSKa0(omega, r, beta, N_SLP, N_lattice)
    
    % Set Quadrature parameters
    num_phi    = 11;    % Odd number for Simpsons rule (N > 10).
    num_points = 31;    % Trapezoidal rule use any N > 30.
    
    % Initialize SLP
    SLP_size = 2 * N_SLP + 1;
    SLP_matrix = zeros(SLP_size , SLP_size);  
    m_range = - N_SLP:N_SLP;       
    n_range = - N_SLP:N_SLP;
    
    % Compute the SLP entries
    for i = 1:SLP_size      
        m = m_range(i);            
        for j = 1:SLP_size
            n = n_range(j);
            term1 =  SLPKummer1(m, n, r, beta, omega);   
            term2 = -SLPKummer2(m, n, r);         
            term3 = -SLPKummer3(m, n, r, num_phi); 
            term4 = -SLPKummer4(m, n, r, num_points); 
            term5 =  SLPKummer5(m, n, omega, N_lattice, beta, r);
            SLP_matrix(i, j) = - (term1 + term2 + term3 + term4 + term5);
            % Remark: minus sign becauce of the different conventions used.
       end
    end
   
end

%% --- Define the Terms of the SLP ---

function result = SLPKummer1(m, n, r, beta, omega)

    C = 1 / (omega^2 + (beta(1)^2 + beta(2)^2));

    if m == 0 && n == 0
        integral_value = C * (2 * pi)^2;
    else
        integral_value = 0; % Orthogonality condition
    end
    
    % Scaling factor from SLP
    result = r / (2 * pi) * integral_value;
end


function result = SLPKummer2(m,n,r)

    C = 1/12 - log(2) / (2 * pi) + 0.5 * r^2;
    
    % Orthigonalit of the integral
    if m == 0 && n == 0
        result = C * 2*pi*r; 

    elseif m == 1 && n == 1
        result =  - 0.5*r^3 *pi;

    elseif n == -1 && m == -1
        result = - 0.5*r^3 *pi;

    else
        result = 0;
    end

end


function  result = SLPKummer3(m, n, r, num_phi)
    if n == m
        if n == 0
            torr = 1e-8;  % absolute error
            N    = 4;     % steps in Gauss method

            % --- Define the integrand ---
                inner_integrand = @(theta, phi) 2*log(sinh(pi * r * (sin(theta) - sin(phi))).^2 + sin(pi * r * (cos(theta) - cos(phi))).^2);

            % --- Vectorised Simpson rule ---

                % Generate phi values
                phi_values = linspace(0, 2*pi, num_phi); 
                h = (phi_values(end) - phi_values(1)) / (num_phi -1);    
                
                % --- Integrating the logarithmic singularity ---
    
                % Evaluate the function in phi such that it is only a function
                % of theta and perform an adapted Gauss integral over theta.
        
                inner_integrals = arrayfun(@(phi) adapgauss(@(theta) inner_integrand(theta, phi), 0, 2*pi, torr, N), phi_values);
                
                % The resulting function in (discrete) phi is smooth continous
                % and periodic (it looks that way...).
    
                % Simpson's rule for the outer integral (over phi)
                weights = 2 * ones(1, num_phi);  % Default weight is 2
                weights(2:2:end-1) = 4;          % Every even index gets weight 4
                weights([1, end])  = 1;          % Endpoints get weight 1
            
                result = h / 3 * sum(weights .* inner_integrals);

            % Apply scaling factor (From the Green's function)
            result = -1 / (8 * pi) * result;
        else
            result = pi / abs(n);
        end
    else
        result = 0;  % Orthogonality condition
    end
   
    % Scaling factor from SLP
    result = r / (2 * pi) * result;

 end


function result = SLPKummer4(m, n, r, num_points)

    if n == 0 && m == 0

        f_x = @(x1, x2) ...
            (cos(2*pi*1*x1)/1) .* (exp(2*pi*1*x2) + exp(-2*pi*1*x2)) ./ (exp(2*pi*1) - 1) + ...
            (cos(2*pi*1*x2)/1) .* (exp(2*pi*1*x1) + exp(-2*pi*1*x1)) ./ (exp(2*pi*1) - 1) + ...
            (cos(2*pi*2*x1)/2) .* (exp(2*pi*2*x2) + exp(-2*pi*2*x2)) ./ (exp(2*pi*2) - 1) + ...
            (cos(2*pi*3*x1)/3) .* (exp(2*pi*3*x2) + exp(-2*pi*3*x2)) ./ (exp(2*pi*3) - 1) + ...
            (cos(2*pi*2*x2)/2) .* (exp(2*pi*2*x1) + exp(-2*pi*2*x1)) ./ (exp(2*pi*2) - 1) + ...
            (cos(2*pi*3*x2)/3) .* (exp(2*pi*3*x1) + exp(-2*pi*3*x1)) ./ (exp(2*pi*3) - 1);% 
            %(cos(2*pi*4*x1)/4) .* (exp(2*pi*4*x2) + exp(-2*pi*4*x2)) ./ (exp(2*pi*4) - 1) + ...
            %(cos(2*pi*4*x2)/4) .* (exp(2*pi*4*x1) + exp(-2*pi*4*x1)) ./ (exp(2*pi*4) - 1);

            % Remark: The trucation error in this sum is O(1e-9).
    
        % --- Define the integrand ---
            integrand = @(theta, phi) f_x(r * (cos(theta) - cos(phi)), r * (sin(theta) - sin(phi))); 
    
        % --- Trapezoidal rule ---

            % Define the grid without double-counting 2*pi
            theta_vals  = linspace(0, 2*pi, num_points + 1);
            theta_vals  = theta_vals(1:end-1);  
            phi_vals    = linspace(0, 2*pi, num_points + 1);
            phi_vals    = phi_vals(1:end-1);     
            
            % Step sizes
            dx = 2 * pi / num_points;  % Step size for theta
            dy = 2 * pi / num_points;  % Step size for phi
            
            % Create a grid
            [Theta, Phi] = meshgrid(theta_vals, phi_vals);
            
            % Evaluate the integrand
            Z = integrand(Theta, Phi);
            
            % Trapezoidal weights
            weights = ones(num_points, 1);
            weights(1)      = 0.5; 
            weights(end)    = 0.5;  
            
            % Compute the integral using proper weighting
            integral_value = sum(sum((weights * weights') .* Z)) * dx * dy;
        
        % --- Apply the scaling ---
            result = (1 / (4 * pi)) * integral_value;
            result = r / (2 * pi) * result;
    else
        result = 0; % orthogonality
    end
end


function lattice_sum_value = SLPKummer5(m, n, w, n_lat, beta, R)

    % --- Initialize the dual lattice vectors ---

        % Dual lattice vectors
        a = [2*pi, 0];
        b = [0, 2*pi]; 
    
        % Precompute terms
        psi_mn = (-1)^n * 1i^(m - n); 
    
        % Create mesh grid for the q values
        [i_vals, j_vals] = meshgrid(-n_lat:n_lat, -n_lat:n_lat);
        i_vals = i_vals(:); 
        j_vals = j_vals(:); 
        
        % Remove the (0,0) element as it is skipped in the sum
        valid_idx = ~(i_vals == 0 & j_vals == 0); 
        i_vals = i_vals(valid_idx);
        j_vals = j_vals(valid_idx);
    
        % Compute the summation points
        q = [i_vals, j_vals] * [a; b];  
    
    % --- Vectorised compuation of the Green's function ---

        % Precompute constants
        beta_norm_sq = beta(1)^2 + beta(2)^2;
        q_norm_sq    = q(:,1).^2 + q(:,2).^2;
    
        % Define the Green's function
        denominator1 = w^2 + beta_norm_sq - 2j * (beta(1) * q(:,1) + beta(2) * q(:,2)) - q_norm_sq;
        denominator2 = q_norm_sq;
        common_denominator = denominator1 .* denominator2;
        green_values = (beta_norm_sq + w^2 - 2j * (beta(1) * q(:,1) + beta(2) * q(:,2))) ./ common_denominator;
    
        % Precompute the angular factor and Bessel functions
        psi_vals = atan2(q(:,2), q(:,1)); 
        scale = 2*pi*R * exp(1i * psi_vals * (n - m)) * psi_mn;  
        arg = R * sqrt(q(:,1).^2 + q(:,2).^2);  
        bessel_values = besselj(m, arg) .* besselj(n, arg); 

    % --- Vectorised lattice sum ---
        res = scale .* bessel_values .* green_values;
        lattice_sum_value = sum(res); 
end


% Non vectorised (i.e. slow lattice sum implmentation)

%{
function green_sum = Green(q, k, beta)
    
    % ----------- Precompute some terms ONLY ONCE -----------
  
    beta_norm_sq = beta(1)^2 + beta(2)^2;  
    q_norm_sq = q(1)^1 + q(2)^2;  

    % ----------- Compute the denominators and numerator -----------
    denominator1 = k^2 + beta_norm_sq - 2j * (beta(1)*q(1) + beta(2)*q(2)) - q_norm_sq;

    denominator2 =  q_norm_sq;

    % ----------- Compute the common denominator -----------
    common_denominator = denominator1 * denominator2;

    % ----------- Compute the Green's function result with common denominator -----------
    green_sum = ( beta_norm_sq + k^2 - 2j * (beta(1)*q(1) + beta(2)*q(2)) ) / common_denominator;  
end

function lattice_sum_value = SLPKummer5(m, n, w, n_lat, beta, R)
    a = [2*pi, 0]; % Initialize the dual lattice
    b = [0, 2*pi];

    lattice_sum_value = 0;
    psi_mn = (-1)^n * 1i^(m - n); % Precompute this term

    for i = -n_lat:n_lat 
        for j = -n_lat:n_lat

            % --- Skip the summation over 0 ---
            if i == 0 && j == 0
                continue; 
            end
            q = i * a + j * b;
            alpha_q = q;
            psi = atan2(alpha_q(2), alpha_q(1));
            scale = 2*pi*R * exp(1i * psi * (n - m)) * psi_mn;
            arg = R*sqrt(alpha_q(1)^2 + alpha_q(2)^2);
            res = scale * besselj(m, arg) * besselj(n, arg) * Green(q, w, beta);
      
            lattice_sum_value = lattice_sum_value + res;
        end
    end    
    lattice_sum_value = 1/(1) * lattice_sum_value; %rescale to volume of the unit cell, (|Y| = 1)
end

%}
