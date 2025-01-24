%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [General Spectral Plot for a 1D resonator Chain]
    ---------------------------------------------------------
%}

% Objective:
%   Generate a spectral plot for a resonator chain in the general regime. 
%   A complete derivation of the formula used in this implementation 
%   can be found in arXiv:2409.14537.

clear all;
close all;

% --- Set the parameters ---
    n = 1.8;                % wave number coupling
    a = 0.2;                % length of resonator
    del = 1;                % contrast | subwavelength regime (del=0.1) 
    L = 1;                  % length of the unit cell
    num_points = 200;       % number of scatter points
    freq_range = 5.5;       % frequency range

% --- Initialise arrays ---
    k_values = linspace(0,freq_range, num_points);
    alpha_results = zeros(length(k_values), 2);
    beta_results  = zeros(length(k_values), 2);

% --- Compute the spectral bands ---
    for i = 1:length(k_values)
        alpha_results(i, :) = alpha(k_values(i), a, del, n, L);
        beta_results(i, :)  = beta(k_values(i), a, del, n, L);
    end

% --- Plotting the spectral Bands ---
    figure;
    lw = 4.5;   % Linewidth
    fs = 18;    % Fontsize
    
    % Plotting of the bands
    pbaspect([8 6 1]); 
    plot(k_values, abs(real(alpha_results(:, 1))) , 'k', 'LineWidth', lw);
    hold on;
    plot(k_values, -abs(real(alpha_results(:, 2))), 'k', 'LineWidth', lw);
    plot(k_values, abs(real(beta_results(:, 1)))  , 'r', 'LineWidth', lw);
    plot(k_values, -abs(real(beta_results(:, 2))) , 'r', 'LineWidth', lw);

    % Adding labels, legend and ticks
    xlabel('k', 'FontSize', fs);
    ylabel('' , 'FontSize', fs);
    legend('\alpha(k)','', '\beta(k)', 'Location', 'southeast', 'FontSize', fs + 2);
    set(gca, 'FontSize', fs + 4); 
    box on;
    saveas(gcf, 'General_Bandplot.pdf', 'pdf');
    hold off;
    

%% --- Generate the Transfer matrix ---

% --- Transfer Matrix over first boundary ---
    function T = T1(k, a, del, n) 
    
        A = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];
        B = [1 + n/del, 1 - n/del; 1 - n/del, 1 + n/del];
        C = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];
    
        T =  0.5 * A * B * C;
    end

% --- Transfer Matrix over second boundary ---
    function T = T2(k, a, del, n) 

        A = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];
        B = [1 + del/n, 1 - del/n; 1 - del/n, 1 + del/n];
        C = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];
    
        T =  0.5 * A * B * C;
    end

% --- Phase shift due to the propagation --- 
    function F = FL(k, L)
    
        F = [exp(-1i * k * L), 0 ; 0, exp(1i * k * L)];
    end 

% --- Transfer matrix over one unit cell --- 
    function M = M(k, a, del, n, L)
    
        M = T1(k, a, del, n) * T2(k, a, del, n) * FL(k,L);
    end

% --- Compute the complex quasimomentum ---
    function b = beta(k, a, del, n, L) 
    
        eigenvalues = eig(M(k, a, del, n, L));
        b = [1 / L * log(abs(eigenvalues(1))), 1 / L * log(abs(eigenvalues(2)))];
    end 

% --- Compute the real quasimomentum ---
    function a = alpha(k, a, del, n, L)
    
        eigenvalues = eig(M(k, a, del, n, L));
        a = [-1/L * angle(eigenvalues(1)), -1/ L * angle(eigenvalues(2))];
    end
