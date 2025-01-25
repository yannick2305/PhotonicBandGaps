%{
    -----------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Evaluate the outside field directly]
    -----------------------------------------------------------
%}

clear all;
close all;

N_multi = 2;        

alpha =  pi * [1,1]; % Keep alpha fixed at [pi, pi]

% --- Precomputed vlaues where SLP has nontrivial kernel ---

%     Remark: other parameters can be precomputed in TwoD_FindKernelSLP.m

%       -> change R below accordingly

    % --- Values for alpha = [pi, pi] ---
    %beta = pi * [-1, 1];            % Dilute resonators R = 0.00001
    %beta = 3.98415  * [1,1];        % R = 0.05  SLP = 2
    %beta = 4.61125 * [1,1];         % R = 0.1   SLP = 2
    beta = 5.49465 * [1,1];         % R = 0.3   SLP = 2 NEW
    %beta = 4.611103 * [1,1];        % R = 0.1   SLP = 4

    % --- Values for alpha = [1, 1] ---
     %beta = 6.13148 * [1, -1];      % R = 0.1; SLP = 4, alpha = [1,1]
     %beta = 1.25195 * [1, -1];      % R = 0.001; SLP = 2, alpha = [1,1]
  

    R =  0.3;  % Resonator radius

% --- Define the parameters ------------------------------------------
    k0 = 0.0001; %0.00001;
    N = 1;
    D = 1;              
    c1 = 1/2*D*[0,0];  
    c = c1;
    N_lattice = 30;       % Use N = 100 due to slow lattice sum convergence
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 0.001; 
    vb = 1;

% Compute the omplex single layer potential 
    matS = makeS(k0,R,alpha,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice);
    matR = makeR(k0, R, alpha, beta, N_multi, N_lattice * 1);
    
    matSR = (-matR + matS)  ;  

   
% --- Comput the pseudokernel of the SLP ---
    K = null(matSR);
    Rank = rank(matSR);
    
    fprintf('\nRank of matSR: %d\n', Rank);
    
    fprintf('\nBasis for the Kernel:\n');
    disp(K);
    
    % Define a tolerance for determining "close to zero"
    tolerance = 0.00001;
    
    singularValues = svd(matSR);
    [U, S, V] = svd(matSR);
    smallSingularIndices = find(singularValues < tolerance);
    
    if ~isempty(smallSingularIndices)
        pseudoKernel = V(:, smallSingularIndices);
        
        % Display the pseudokernel
        disp('Pseudokernel (Approximate Null Space) of SR:');
        disp(pseudoKernel);
    else
        disp('No pseudokernel (all singular values are above tolerance).');
    end
    
    % Double check the computed vector lies in the pseudokernel
    check = (-matR + matS ) * pseudoKernel;


%% --- Compute the outside Field solution ---

ker = pseudoKernel(:,1); % Select one of the pseudokernel vectors

% --- Plot the outside field  ---
    % Define grid values
    x1_vals = linspace(-0.5, 1.5, 50);
    x2_vals = linspace(-0.5, 1.5, 50);
    [x1, x2] = meshgrid(x1_vals, x2_vals);

    % Initialize storage for results
    Res = zeros(size(x1));

    % Loop over the grid
    for i = 1:size(x1, 1)
        for j = 1:size(x2, 2)
            x = [x1(i, j), x2(i, j)];
            
            Res(i, j) = evaldFieldR(x, ker, alpha, beta, k0, R, N_multi, N_lattice);
        end
    end

 
%% --- Plot the outside field ---

Sol = real(Res);
     
lim = 0.1;  % Size of the plotting window
fs  = 26;   % Fontsize of the annotations

% --- Create a surface plot ---
    figure;
    surf(x1, x2, (Sol));
    shading interp;
    colorbar;
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fs + 4); 
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fs + 4); 
    set(gca, 'FontSize', fs-4); 
    %axis equal;
    axis tight;
    view(2);
 
    clim([-lim, lim]); 
    zlim([-lim, lim]); 


%% --- Define the functions ---

function sol = evaldFieldR(x, ker, alpha, beta, w, R, N_SLP, N_lattice)

    % Number of elements for SLP_vec
    SLP_size = 2 * N_SLP + 1;
    n_range  = -N_SLP:N_SLP;

    % Initialize SLP_vec
    SLP_vec = zeros(1, SLP_size);

    % Loop over the range of n and compute the LatticeSum for each n
    for i = 1:SLP_size
        n = n_range(i);
        SLP_vec(i) = LatticeSum(n, x, alpha, beta, w, R, N_lattice);
    end

    % Compute the solution as the dot product
    sol =  SLP_vec * ker;
end


function res = LatticeSum(n, x, alpha, beta, w, R, N_lattice)
    % Generate lattice indices
    [a, b] = meshgrid(-N_lattice:N_lattice, -N_lattice:N_lattice);
    q = 2 * pi * [a(:), b(:)]; % Each row of q is a vector in the reciprocal lattice
    
    % Compute q_alpha for all lattice points
    q_alpha = q + alpha; % Each row corresponds to alpha + q
    
    % Compute norms and dot products
    q_alpha_norm     = sqrt(sum(q_alpha.^2, 2));    % Vector magnitude |alpha + q|
    q_alpha_dot_x    = q_alpha * x';                % Dot product (alpha + q) · x
    q_alpha_dot_beta = q_alpha * beta';             % Dot product (alpha + q) · beta
    q_alpha_norm_sq  = q_alpha_norm.^2;             % Squared norm |alpha + q|^2
    beta_norm_sq     = sum(beta.^2);                % Norm of beta squared

    arg_q_alpha = atan2(q_alpha(:, 2), q_alpha(:, 1)); 

    % Compute the Green function for all lattice points
    denominator1 = w^2 + beta_norm_sq  - 2j * q_alpha_dot_beta - q_alpha_norm_sq ;
    green_sum    = 1 ./ (denominator1);
    
    % Compute the Bessel function arguments and values
    argument      = R * q_alpha_norm;       % Argument for besselj
    bessel_values = besselj(n, argument);   % Compute Bessel function
    
    % Combine terms to compute the lattice sum
    exponential_term = exp(1i * q_alpha_dot_x); % e^(i (alpha + q) · x)
    lattice_sum      = sum(exponential_term .* bessel_values .* exp(1i * n * arg_q_alpha)  .* green_sum); % Summation over all q
    
    % Scale the result
    scale = 2 * pi * R * (-1i)^(n);
    res   = scale * lattice_sum;
end


