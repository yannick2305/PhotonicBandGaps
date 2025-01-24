%{
    ------------------------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Generate complex single layer potential  without accelerated lattice sums]
    ------------------------------------------------------------------------------------------
%}

function SLP_matrix = makeRIntuitive(w, R, alpha, beta, N_SLP, N_lattice)
    %
    % Description:
    %   This function computes the remainder Single Layer Potential (RSLP) 
    %   matrix by performing the boundary integral.
    %
    % Input Arguments:
    %   w          - Frequency (scalar).
    %   R          - Resonator radius (scalar).
    %   alpha      - Real part of the quasimomentum (R^2-vector).
    %   beta       - Complex part of the quasimomentum, (R^2-vector).
    %   N_SLP      - Size of the Single Layer Potential (SLP) matrix (scalar).
    %   N_lattice  - Number of lattice points used in the Green's function (scalar).
    %
    % Output Arguments:
    %   SLP_matrix - The resulting Single Layer Potential matrix (size: 2N+1).

    
    % ----------- Precompute some terms ONLY ONCE -----------
        % Initialize SLP
        SLP_size = 2 * N_SLP + 1;
        SLP_matrix = zeros(SLP_size , SLP_size);  
        [i_vals, j_vals] = meshgrid(-N_lattice:N_lattice, -N_lattice:N_lattice);  
    
        % Dual lattice vectors
        a = [2*pi, 0]; 
        b = [0, 2*pi];
        
        % Compute the full q grid (precompute q1, q2, and q)
        q1 = i_vals * a(1);  % Dual lattice
        q2 = j_vals * b(2);  % Dual lattice
        q = [q1(:), q2(:)];  % Reshape q to a 2D array (2N x 2)
        
        % Define the range of m and n
        m_range = -N_SLP:N_SLP;   % Finite Fourier Basis
        n_range = -N_SLP:N_SLP;
        k =  -2*N_SLP:2*N_SLP;

        % Precompute the terms involed in the sum for the Green's function
        green = Green(q, w, alpha, beta); 
    
        % Precompute the argument for the Bessel functions used in the latice sum
        alpha_q = alpha + q;  
        arg = R * sqrt(alpha_q(:,1).^2 + alpha_q(:,2).^2);  % Argument for Bessel
        [unique_arg, ~, idx_map] = unique(arg);             % Precompute unique arguments and indices
        
        % Precompute psi, which only depends on alpha_q
        psi = atan2(alpha_q(:,2), alpha_q(:,1)); 
        
        % Precompute the circumference
        TWOpiR = 2 * pi * R;

        % Precompute for lattice sum
        exp_result = TWOpiR * exp(1i * psi * k)   .* green; % .* (1i.^k)
        
        % Precompute Bessel functions
        max_order = max(abs([m_range, n_range])); 
        precomputed_bessel = cell(max_order + 1, 1);
        
        for order = 0:max_order
            precomputed_bessel{order + 1} = besselj(order, unique_arg);
        end
        
    % ----------- Pass precomputed terms into Lattice sum -----------
        % Evaluate entries of the single layer potential
        for i = 1:SLP_size          
            m = m_range(i);        
            for j = 1:SLP_size
                n = n_range(j);
                % Precompute lattice sum arguments
                kind = m - n + (2*N_SLP + 1);
                if mod(n, 2) == 0  
                    scale =  exp_result(:, kind );   % Assign +1 for even n
                else
                    scale =  -exp_result(:, kind );   % Assign -1 for odd n
                end
                % Perform the lattice sum
                SLP_matrix(i, j) = (1i)^(m+n) * lattice_sum(m, n, idx_map, scale, precomputed_bessel);
            end
        end

end


function green_sum = Green(q, w, alpha, beta)
    %
    % Description:
    %   This function computes the Single Layer Potential (SLP) matrix by
    %   performing the boundary integral.
    %
    % Input Arguments:
    %   w          - Frequency (scalar).
    %   q          - Lattice vector (R^2-vector).
    %   alpha      - Real part of the quasimomentum (R^2-vector).
    %
    % Output Arguments:
    %   green_sum  - Term involved in lattice sum of Green's function (scalar).


    % ----------- Precompute some terms ONLY ONCE -----------
        q_alpha          = alpha + q;  
        q_alpha_norm_sq  = sum(q_alpha.^2, 2); 
        q_alpha_dot_beta = sum(q_alpha .* beta, 2); 
        beta_norm_sq     = sum(beta.^2);  

    % ----------- Compute the denominators and numerator -----------
        denominator1 = q_alpha_norm_sq + 2j * q_alpha_dot_beta - w^2 - beta_norm_sq;
        %denominator2 = q_alpha_norm_sq - w^2;
        %numerator = beta_norm_sq - 2j * q_alpha_dot_beta;

    % ----------- Compute the Green's function result -----------
        green_sum = 1 ./ (denominator1);  
end


function lattice_sum_value = lattice_sum(m, n, idx_map, scale, precomputed_bessel)
    %
    % Description:
    %   This function computes the Single Layer Potential (SLP) matrix by
    %   performing the boundary integral.
    %
    % Input Arguments:
    %   m                   - Row index of Single Layer Potential (scalar).
    %   n                   - Column index of Single Layer Potential (scalar).
    %   green               - Precomputed term in lattice sum (scalar).
    %   scale               - Precomputed scale (vector).
    %   idx_map             - Precomputed argument index (vector).
    %   precomputed_bessel  - Precomputed bessel functions (vector).
    %
    % Output Arguments:
    %   lattice_sum_value  - Lattice sum defining the rmainder SLP (scalar).


    % ----------- Evaluate precomputed Bessel function -----------
        if m >= 0
            bessel_m_unique = precomputed_bessel{m + 1};
        else
            % Symmetry: J_{-m}(x) = (-1)^m J_m(x)
            bessel_m_unique = (-1)^m * precomputed_bessel{-m + 1};
        end
    
        if n >= 0
            bessel_n_unique = precomputed_bessel{n + 1};
        else
            % Symmetry: J_{-n}(x) = (-1)^n J_n(x)
            bessel_n_unique = (-1)^n * precomputed_bessel{-n + 1};
        end
    
        % Map unique Bessel values to the specific indices
        bessel_m = bessel_m_unique(idx_map);
        bessel_n = bessel_n_unique(idx_map);

    % ----------- Perform the lattice sum -----------
        res = scale .* bessel_m .* bessel_n ;  % Element-wise product
        lattice_sum_value = sum(res(:));       % Sum over all elements
end
