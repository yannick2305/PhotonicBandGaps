%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Illustrate the alpha shift within the bandgap]
    --------------------------------------------------------------
%}


%% ------- Create Capacitance matrix ONCE at the beginning -------
%
% Objective: 
%   Compute the capacitance matrix for a 2D array of resonators. 
%   Since a finite capacitance matrix cannot be defined in a truly infinite 
%   2D system, we compute a truncated version of the inverse Floquet-transformed 
%   quasiperiodic capacitance. This truncated matrix serves as a reasonable 
%   approximation for analyzing and measuring the decay behavior in the system.
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
    fprintf('----------------------------\n');
    fprintf('Computation completed in %.4f seconds.\n', elapsed);

%% ------- Compute the decay length and quasimomentum alpha -------
%
% Objective: 
%   This function computes the complex quasimomentum (beta) and the 
%   complex Floquet transform for a given frequency (w) within the bandgap.
%   Additionally, it generates a plot of the density of the real quasimomentum.
%
% Input: 
%   - w: The frequency (w) within the bandgap
%
% Output: 
%   - beta: The complex quasimomentum corresponding to the input frequency.
%   - FloquetTransform: The complex Floquet transform related to the system's 
%                       periodicity and the input frequency.
%   - A plot showing the density of the real part of the quasimomentum.


% --- Set the frequency range within the Bandgap ---
    w = 1.0;        % Bandgap starts at: w0 = 0.579973. Use w > w0.
    Nalpha = 100;   % Discretization of alpha grid. Use Nalpha = 100.

% --- Compute the decay and quasimomentum for gap frequency w ---

    % --- Initialise vectors ---
    I = eye(n^2);
    E = C - w^2 * I;

    % ---  Invert the matrix E ---
    matrix_inverse = E \ I;

    % --- Extract the middle column ---
    i =  floor(n^2/2) + 1;
    ith_column_vector = matrix_inverse(:, i);
    
    % --- Get the defect mode u(x) at matrix coordinates x = (i, j) ---
    indices = repmat(1:n, 1, n);                % Repeat the indices 1 to n, n times
    values = reshape(ith_column_vector, n, n);  % Reshape the vector to an n x n matrix

    % --- Get the decay via the slope of the surface plot ---
    decay_beta = real(log(values(floor(n/2)+1, 1))-log(values(floor(n/2)+1, floor(n/2)))) / (floor(n/2)-1);
    
    % --- Debugging for adjusting the slope values ---
    %plotMatrixSurface(values) % Plot the defect mode in a log scalae.

    % --- Compute the quasimomentum alpha ---
    x = linspace(-pi, pi, Nalpha);         
    y = linspace(-pi, pi, Nalpha);          
    [X, Y] = meshgrid(x, y);                
    functionValues = zeros(Nalpha, Nalpha); 
    
    % --- Evaluate the complex Floquet transform for alpha in grid ---
    parfor i = 1:Nalpha
        for j = 1:Nalpha
            sumValue = 0;
            % Current lattice coordinate r = [X(i, j), Y(i, j)]
            r = [X(i, j)+ 1i * decay_beta, Y(i, j) + 0]; %alpha + i beta
            % Loop over all wavevectors
            for inx = 1:n
                for iny = 1:n
                    % Compute the dot product between r and wavevector k
                    dotProduct = r * shiftCoordinates(inx, iny, n)'; % r dot [a_k, b_k]
                    % Add the term to the sum
                    sumValue = sumValue + values(inx, iny) * exp(1i * dotProduct);
                end
            end
            % Store the result in the function values matrix
            functionValues(i, j) = sumValue;
        end
    end
    
    % --- Determine the peaks of the surface plot (highest alpha density) ---
    [maxValue, maxIndex] = max(functionValues(:));          % Find maximum value and its index
    [row, col] = ind2sub(size(functionValues), maxIndex);   % Convert index to row and column
    maxCoordinates = [x(row), y(col)];                      % Get the coordinates of the maximum alpha value

    % --- Display the numerical Values for alpha (i.e. peaks of surfaceplot) ---
    fprintf('----------------------------\n');
    fprintf('Frequency:   \x03C9 = %.4f\n', w);  
    fprintf('Decay:       \x03B2 = (%.4f, 0)\n', decay_beta);
    fprintf('Quasimoment: \x03B1 = (%.4f, %.4f)\n', maxCoordinates(1), maxCoordinates(2));

    % --- Plotting the alpha density  ---
    figure;
    surf(x, y, abs(functionValues')); 
    xlabel('$\alpha_1$', 'Interpreter', 'latex', 'FontSize', 28);
    ylabel('$\alpha_2$', 'Interpreter', 'latex', 'FontSize', 28);
    zlabel('$|| \hat{u} ||_{L^2}$', 'Interpreter', 'latex', 'FontSize', 28);
    view(3);

    xticks([0, pi]); 
    yticks([0, pi]); 
    xticklabels({'0', '$\pi$'}); 
    yticklabels({'0', '$\pi$'}); 
    set(gca, 'FontSize', 25);
    set(gca,'TickLabelInterpreter','latex')
    
    % --- Save the plot as a PDF ---
    filename = sprintf('alpha_freq_%.2f.pdf', w);   % Add frequency into filename
    %exportgraphics(gcf, filename);                  % Uncomment to Save as PDF file.
  

%% --- Create Movie for the alpha shift --------------------
%
% Input: 
%   - w_values: A range of frequencies within the bandgap.
%
% Output: 
%   - Movie: A visualisation (movie) that illustrates how the real quasimomentum 
%     evolves over the Brillouin zone as the frequency is varied within the 
%     bandgap. The movie shows the dynamic change in quasimomentum as the frequency 
%     moves through the specified range.


% --- Set the frequency range within the Bandgap ---
    w_values =  0.579973 : 0.01 : 2.2;  % Bandgap starts at: w0 = 0.579973;
    Nalpha = 40;                        % Discretization of alpha grid. Use Nalpha = 50.

% --- Create VideoWriter to save the movie ---
    v = VideoWriter('Alpha_Shift_surface_plot.mp4', 'MPEG-4');
    v.FrameRate = 8;    % Set frame rate (Defalt: 10 fps)
    open(v);            % Open video for writing

% --- Initialise vectors ---
    numIterations = length(w_values);
    combinedResults = zeros(numIterations, 2);
    I = eye(n^2);
    figure;

% --- Compute the Complex Floquet transform for frequency range ---
    for idx = 1:numIterations
        w = w_values(idx);
        E = C - w^2 * I;
        matrix_inverse = E \ I;
    
        % --- Extract the middle column of the inverse matrix ---
        i = floor(n^2/2) + 1;
        ith_column_vector = matrix_inverse(:, i);   % Extract middle column
        indices = repmat(1:n, 1, n);                % Repeat the indices 1 to n, n times
        values = reshape(ith_column_vector, n, n);  % Reshape to n x n matrix
    
        % --- Compute the decay length ---
        result = real(log(values(floor(n/2)+1, 1)) - log(values(floor(n/2)+1, floor(n/2)))) / (floor(n/2)-1);
        combinedResults(idx, :) = [result, w];
    
        % --- Quasimomentum computation ---
        x = linspace(0, pi, Nalpha);    % Adjust range from -pi to pi to plot entire Brilloin zone
        y = linspace(0, pi, Nalpha);
        [X, Y] = meshgrid(x, y);
        
        % --- Evaluate the complex Floquet transform for alpha in grid --- 
        r1 = (X(:) + 1i * result);  % (Nalpha^2, 1)
        r2 = Y(:);                  % (Nalpha^2, 1)
        
        [kx, ky] = meshgrid(1:n, 1:n); 
        kVectors = shiftCoordinates(kx, ky, n); 
        
        kVectors = reshape(kVectors, [], 2);
        
        valuesVec = reshape(values, [], 1);
        
        dotProducts = r1 * kVectors(:,1).' + r2 * kVectors(:,2).'; 
        
        expTerms = exp(1i * dotProducts); 
        
        functionValuesVec = expTerms * valuesVec; 
        
        % Reshape back to (Nalpha, Nalpha)
        functionValues = reshape(functionValuesVec, Nalpha, Nalpha);
        
        % --- Surface plot (use abs for visualization) ---
        surf(x, y, abs(functionValues'), 'EdgeColor', 'none');
        xlabel('$\alpha_1$', 'Interpreter', 'latex', 'FontSize', 28);
        ylabel('$\alpha_2$', 'Interpreter', 'latex', 'FontSize', 28);
        zlabel('$|| \hat{u} ||_{L^2}$', 'Interpreter', 'latex', 'FontSize', 28);
        title(['Frequency: $\omega = ', num2str(w, '%.3f'), '$'], 'Interpreter', 'latex');
        view(2);    %Top view
        xticks([0, pi]);
        yticks([0, pi]);
        xticklabels({'0', '$\pi$'});
        yticklabels({'0', '$\pi$'});
        set(gca, 'FontSize', 25);
        set(gca, 'TickLabelInterpreter', 'latex');
        axis equal;
        axis tight;
    
        % --- Capture the frames and write it to the video ---
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    close(v); % Stop the recording


%% --- Defining functions -------------------------------------

function integral_value = Cm(i, j)

    % --- Set values for numerical computation ----------------
    Nq          = 20;   % Convergence   O(dx^2 + dy^2), use 15.
    N_multipole = 1;    % 1 for single resonator (5 for dimer).
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

    JHdata   = makeJHdata0(k0,R,N_multipole);
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
    parfor j = 0:n-1
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
