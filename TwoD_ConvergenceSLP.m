%{
    --------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Convergence rates lattice sums single layer potntial.]
    ---------------------------------------------------------------------------
%}

clear all;
close all;

% --- Set the Parameters ---
    R = 0.05;       % Resonator radius.
    N_SLP = 2;      % Dimension of SLP (use N_SLP = 2).
    k0 = 0.00001;   % Frequency.

    % Generic values for alpha and beta
    beta  = [1,   0]; 
    alpha = [pi, pi];

% --- Generate the "true" reference SLP using the three methods  ---
    Lattice_ref = 1000; 
    Kummer    =  makeSKa0(k0, R, beta, N_SLP, Lattice_ref);
    Intuitive =  makeRIntuitive(k0, R, alpha, beta, N_SLP, Lattice_ref);
    Old       = -makeR(k0, R, alpha, beta, N_SLP, Lattice_ref); 

% --- Initialise arrays ---
    n_values = round(logspace(log10(1), log10(500), 30)); 
    convergence1 = zeros(length(n_values), 1);
    convergence2 = zeros(length(n_values), 1);
    convergence3 = zeros(length(n_values), 1);

% --- Loop over the truncation parameter n ---
    parfor idx = 1:length(n_values)
        n = n_values(idx);
        
        A1_n =  makeSKa0(k0, R, beta, N_SLP, n);
        A2_n =  makeRIntuitive(k0, R, alpha, beta, N_SLP, n);
        A3_n = -makeR(k0, R, alpha, beta, N_SLP, n);
        
        % Difference between reference matrix and its approximation
        convergence1(idx) = norm(Kummer    - A1_n, 'fro');
        convergence2(idx) = norm(Intuitive - A2_n, 'fro');
        convergence3(idx) = norm(Old       - A3_n, 'fro');
    end

%% --- Compute and plot the decay rate ---

% --- Compute the algebraic decay rate ---
    log_n = log(n_values);  
    log_convergence1 = log(convergence1);
    log_convergence2 = log(convergence2);
    log_convergence3 = log(convergence3);
    
    % Fit a linear model to log(n) vs log(convergence)
    p1 = polyfit(log_n(2:end), log_convergence1(2:end), 1); 
    p2 = polyfit(log_n(2:end), log_convergence2(2:end), 1);
    p3 = polyfit(log_n(2:end), log_convergence3(2:end), 1);
    
    % Extract the exponents from the slope of the fit line (slope = -p)
    p1_value = -p1(1); 
    p2_value = -p2(1);
    p3_value = -p3(1);
    
    % Display the result
    fprintf('-------------------------------------\n');
    disp(['Decay exponent for Kummer:  ', num2str(p1_value)]);
    disp(['Decay exponent for Direct:  ', num2str(p2_value)]);
    disp(['Decay exponent for Classic: ', num2str(p3_value)]);

% --- Plot the truncation error ---
    figure('Position', [100, 100, 1000, 400]); 
    loglog(n_values, convergence1, 'LineWidth', 4, 'LineStyle', '-'  , 'DisplayName', 'Kummer' );
    hold on;
    loglog(n_values, convergence2, 'LineWidth', 4, 'LineStyle', '- -', 'DisplayName', 'Direct' );
    loglog(n_values, convergence3, 'LineWidth', 4, 'LineStyle', ':'  , 'DisplayName', 'Classic');
    hold off;
    
    xlabel('$n$ (Truncation Parameter)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Absolute error', 'Interpreter', 'latex', 'FontSize', 16);
    legend('Location', 'northeast', 'FontSize', 14, 'Interpreter', 'latex');
    
    set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
    grid on;
    
    % Save the figure as a PDF
    %print(gcf, 'ConvergenceRateSLP.pdf', '-dpdf', '-bestfit');
