%{
    -------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Compute the decay length for bandgap resonant frequency]
    -------------------------------------------------------------------------
%}

%% ------- Create Capacitance matrix ONCE at the beginning -------
%
% Objective: 
%   Compute the capacitance matrix for a 2D array of resonators. 
%   Since a finite capacitance matrix cannot be defined in a truly infinite 
%   2D system, we compute a truncated version of the inverse Floquet-transformed 
%   quasiperiodic capacitance. This truncated matrix serves as a reasonable 
%   approximation for analyzing and measuring the decay behaviour in the system.
%
% Remark 1: 
%   It is recommended to utilize MATLAB's Parallel Computing Toolbox 
%   to enhance performance and reduce computation time.

clear all;
close all;
tic; 

% --- Set the size of the capacitance matrix ---
    n = 11;  % nxn system of resonators n odd.
             % n = 7, 9, 11, 13, 15, 17  produce good results.
             % Runtime O(n^2).

% Remark 2: 
%   Since we are analyzing exponential decay, the computed values may 
%   approach MATLAB's machine epsilon (eps = 2.2204e-16). To maintain 
%   numerical accuracy and avoid issues with extremely small values, 
%   we restrict the calculations to values greater than eps, which 
%   corresponds to considering n < 18.

% --- Generate the truncated capacitance matrix --- 
    C = generateBlockToeplitzMatrix(n);
    C = real(C); % Up to numerical error the Capacitance entries are real.

% --- Measure the Computation time ---
    elapsed = toc;  
    fprintf('------------------------------------------\n');
    fprintf('Computation completed in %.4f seconds.\n', elapsed);

%% ---- Compute the decay length for bandgap frequencies ---

% --- Define the system parameters---
    % Set the frequency range within the Bandgap 
    w_values =  0.579973 : 0.02 : 1.5;  % Bandgap starts at: w0 = 0.579973

% --- Initialize variables ---
    numIterations = length(w_values);
    combinedResults = zeros(numIterations, 2);
    I = eye(n^2);

% --- Iterate over the frequencies  ---
    for idx = 1:numIterations
        % --- Initialise vectors ---
        w = w_values(idx);
        E = C - w^2 * I;
        
        % ---  Invert the matrix E ---
        matrix_inverse = E \ I;
    
        % Extract the middle column of the inverse matrix, others are same but shifted
        i =  floor(n^2/2) + 1;
        ith_column_vector = matrix_inverse(:, i);
        
        % --- Create an array of indices ---
        indices = repmat(1:n, 1, n);                % Repeat the indices 1 to n, n times
        values = reshape(ith_column_vector, n, n);  % Reshape the vector to an n x n matrix
        
         % --- Get the decay via the slope of the surface plot ---
        result = real(log(values(floor(n/2)+1, 1))-log(values(floor(n/2)+1, floor(n/2)))) / (floor(n/2)-1);
        combinedResults(idx, :) = [result, w];
    
    end

% --- Plot the decay as a function of the frequency ---
    figure; 
    plot(combinedResults(:,1), combinedResults(:,2), 'x', 'MarkerSize', 5, 'LineWidth', 2, 'Color', 'b');
    xlabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 28);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 28);
    set(gca, 'FontSize', 20);

% --- Display the results ---
    %{
    fprintf('------------------------------------------\n');
    disp('  [Decay,   Frequency]');
    for i = 1:size(combinedResults, 1)
        fprintf('  %.5f   %.3f\n', combinedResults(i, 1), combinedResults(i, 2));
    end
    %}


%% --- Defining functions -------------------------------------

function integral_value = Cm(i, j)

    % --- Set values for numerical computation ----------------
    Nq          = 20;   % Convergence   O(dx^2 + dy^2), use 15.
    N_multipole = 1;    % 1 for single resonator, 5 for dimer.
    N_lattice   = 5;    % lattice sum is exponentially convergent use 3-4.

    % --- Define the Brillouin zone ---------------------------
    a = -pi;
    b =  pi;

    % --- Parameters of the resonator lattice ---------
    k0  = 0.00001;          
    R   = 0.05; 
    vol = pi*R^2;
    d_zeta=makezetadata;
    D = 1; 
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x,L2y];
    
    % --- Single resonator parameters ---
    c1 = 1/2*D*[1,1];
    N = 1; c = [c1];

    % --- Dimer resonator parameters ---
    %c1 = D*[1/2,1/2];
    %c2 = D*[3/4,1/2];
    %N = 2; c = [c1;c2];

    JHdata = makeJHdata0(k0,R,N_multipole);
    JHijdata = makeJHijexpdata(k0,c,N_multipole);

    % --- Material parameters for bubble ---
    rho_0 = 1000;
    rho_b = 1;
    kappa_0 = 1000;
    kappa_b = 1;
    v0 = sqrt(kappa_0/rho_0);
    vb = sqrt(kappa_b/rho_b);
    delta = 1e-3; % =rho_b/rho_0;

    % --- Implmentation of the Trapezoidal rule --------------
    f = @(x, y) ((delta * vb^2)/vol) * makeC(k0, R, [x, y], L1x, L2, d_zeta, JHdata, JHijdata, N, N_multipole, N_lattice) * exp(-1i*(x*i + y*j));
 
    % Create x and y vectors with n points from a to b
    x = linspace(a, b, Nq);
    y = linspace(a, b, Nq);
    
    % Size of the grid
    nx = length(x);
    ny = length(y);
    
    % Create meshgrid for function evaluation
    [X, Y] = meshgrid(x, y);
    
    % Evaluate the function at each point in the grid
    Z = arrayfun(f, X, Y);
    
    % Get the spacing between points in x and y directions
    dx = (b - a) / (Nq - 1);
    dy = (b - a) / (Nq - 1);
    
    % Compute the integral using the trapezoidal rule
    % Add weights for the trapezoidal rule at the boundaries
    weights_x = ones(1, nx);
    weights_y = ones(1, ny);
    
    % Half weight for boundary points
    weights_x(1)    = 0.5;
    weights_x(end)  = 0.5;
    weights_y(1)    = 0.5;
    weights_y(end)  = 0.5;
    
    % Apply trapezoidal rule to sum over x and y (vectorised sum)
    integral_value = (1/(2*pi)^2) * dx * dy * sum(sum(Z .* (weights_y' * weights_x)));
end

function T = generateToeplitz(i, n)
    % Generate the first row and the first column for the Toeplitz matrix
    first_row = zeros(1, n);
    first_col = zeros(n, 1);
    
    % Populate the first row and first column using the function C(i, j)
    for j = 0:n-1
        first_row(j+1) = Cm(i, j);
        first_col(j+1) = Cm(i, j);
    end
    
    % Create the Toeplitz matrix using the first row and column
    T = toeplitz(first_col, first_row);
end

function B = generateBlockToeplitzMatrix(n)
    % Precompute a vector of blocks
    T_blocks = cell(1, n);
    parfor shift = 0:n-1
        T_blocks{shift+1} = generateToeplitz(shift, n);
    end

    % Create a block Toeplitz matrix
    % Start with a scalar Toeplitz pattern
    scalar_toeplitz = toeplitz(0:n-1);

    % Replace each scalar with the corresponding block
    B = cell2mat(arrayfun(@(s) T_blocks{s+1}, scalar_toeplitz, ...
                          'UniformOutput', false));
end

function plotMatrixSurface(A)
    % A: Input square matrix (n x n)
    A = real(log(A));
    % Get the size of the matrix
    [nRows, nCols] = size(A);
    
    % Create a grid for the X and Y coordinates from 0 to n-1
    [X, Y] = meshgrid(0:nCols-1, 0:nRows-1);  % Coordinates from (0,0) to (n-1,n-1)
    
    % Create a surface plot
    figure; % Create a new figure window
    surf(X, Y, A, 'EdgeColor', 'interp'); % Use 'interp' for smooth edges
    
    % Set labels and title
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 20);
    zlabel('$\log(|u(x)|)$', 'Interpreter', 'latex', 'FontSize', 20);

    %title('Surface Plot of Matrix Entries');
    
    % Adjust view angle for better visualization
    view(0, 90);  % Top view to emphasize grid structure
    
    % Optional: Add color bar
    colorbar;
    
    % Set the ticks to be integers corresponding to matrix dimensions
    set(gca, 'XTick', 0:nCols-1, 'YTick', 0:nRows-1);
    
    % Set the axis limits for clarity
    xlim([-0.5, nCols-0.5]);
    ylim([-0.5, nRows-0.5]);
    zlim([min(A(:)), max(A(:))]);
    
    % Enhance grid visibility
    grid on;
end

function shiftedVector = shiftCoordinates(i, j, n)
    % Function to shift the normal matrix coordinates (i, j) to centered coordinates
    % Inputs:
    %   i - Row index in range [1, n]
    %   j - Column index in range [1, n]
    %   n - Size of the matrix (assumed to be square n x n)
    
    % Compute the center of the matrix
    center = ceil(n / 2);
    
    % Shift the coordinates and return as a vector [a, b]
    shiftedVector = [i - center, j - center];
end
