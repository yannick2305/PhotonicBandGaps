%{
    -------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Compute the field solution, with 0 multipole.]
    -------------------------------------------------------------
%}

clear all;
close all;

% --- Fixed Parameters --- 
    a0    = 1;          % SLP coeficient
    Y     = 1;          % Volume of first Brillouin zone
    k     = 0.01;       % Frequency
    alpha = [pi, pi];   % Quasimomenta, has to be fixed to M
    N = 60;            % Lattice summation size
    Resolution = 100;   % Resolution of the solution u(x)
    size = 1.5;         % Size of the surface plot [0.5, size] x [0.5, size]

% --- Variable Parameters --- 
    r = 0.05;            % Resonator radius

    beta = 4.09699 * [1, 1];   % set R = 0.05
    %beta = 6.88444 * [1, 1];   % set R = 0.3

% Remark: 
%       β is precomputed such that the SLP has non-trivial kernel, that
%       way the solution on the boundary of the resonator is equal to 0.


% --- Compute the lattice sum ---

    % --- Precompute what is independent of x ---

    % Generate the dual lattice
    [q_x, q_y] = meshgrid((-N:N) * 2*pi, (-N:N) * 2*pi); % Dual lattice
    q_vectors = [q_x(:)'; q_y(:)']; % Combine into a 2x(2N+1)^2 array
    
    % Precompute alpha + q and norms
    alpha_q = q_vectors + alpha(:); % 2xK array
    norm_alpha_q = sqrt(sum(alpha_q.^2, 1)); % Norm for each column
    
    % Compute dot(beta, alpha + q) for the denominator
    dot_beta_alpha_q = beta(:)' * alpha_q; % Row vector of dot products
    
    % Denominator 
    denom = k^2 + norm(beta)^2 - 2i * dot_beta_alpha_q - norm_alpha_q.^2;
    
    % Bessel function term 
    bessel_term = besselj(0, r * norm_alpha_q);
    
    % Prefactor 
    prefactor = (2 * pi * r * a0 / abs(Y));

% --- Find solution at Grid points ---

    x1_range = linspace(0.5, size, Resolution); 
    x2_range = linspace(0.5, size, Resolution); 
    [x1, x2] = meshgrid(x1_range, x2_range); 
    
    % Initialize result matrix
    u_grid = zeros(Resolution);
    
    % Evaluate u(x) at each point on the x grid
    parfor idx = 1:numel(x1)

        x_eval = [x1(idx); x2(idx)];
        
        dot_alpha_q_x = sum(alpha_q .* x_eval, 1); % 1xK dot products
        
        exp_term = exp(1i * dot_alpha_q_x);
        
        summand = exp_term .* bessel_term ./ denom; % Vectorised summation
        u_grid(idx) = sum(summand); 
    end
    
    u_grid = prefactor * u_grid;


%% --- Plot the outside solution ---

fs = 26;    % Define the label fontsize

% --- Plotting ---
    figure;
    surf(x1, x2, real(u_grid), 'EdgeColor', 'none');
    colorbar;
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fs + 4); 
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fs + 4); 
    set(gca, 'FontSize', fs-4); 
    axis equal;
    axis tight;
    view(2);

%% --- Check correctness of the outside solution ---

U = real(u_grid);   % real part of the solution (imaginary part should be 0)

% Plug the solution U into the original PDE
PDE =  PDEoperator(U, x1_range, x2_range, beta, k);

fs = 26;    % Define the label fontsize

% --- Plotting ---
    figure;
    surf(x1_range, x2_range, abs(PDE));
    shading interp;
    colorbar;
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fs + 4); 
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fs + 4); 
    set(gca, 'FontSize', fs-4); 
    axis equal;
    axis tight;
    view(2);
    
    clim([0, 0.2]);
    zlim([-4, 4]); 
    view(2);

%% --- Differential operator using FEM ---

function result = PDEoperator(v, x, y, beta, k)

    % Constants
    beta_x = beta(1);
    beta_y = beta(2);
    beta_norm_sq = beta_x^2 + beta_y^2;
    k_sq = k^2;

    % Uniform grid spacing
    dx = x(2) - x(1);
    dy = y(2) - y(1);

    % Laplacian: ∆v = d²v/dx² + d²v/dy²
    d2v_dx2 = (circshift(v, [-1, 0]) - 2 * v + circshift(v, [1, 0])) / dx^2;
    d2v_dy2 = (circshift(v, [0, -1]) - 2 * v + circshift(v, [0, 1])) / dy^2;
    laplacian = d2v_dx2 + d2v_dy2;

    % Gradient: ∇v = [dv/dx, dv/dy]
    dv_dx = (circshift(v, [-1, 0]) - circshift(v, [1, 0])) / (2 * dx);
    dv_dy = (circshift(v, [0, -1]) - circshift(v, [0, 1])) / (2 * dy);

    % Dot product: 2 * β · ∇v
    beta_dot_grad_v = 2 * (beta_x * dv_dx + beta_y * dv_dy);

    % Final operator
    result = laplacian - beta_dot_grad_v + (k_sq + beta_norm_sq) * v;
end

