%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Subwavelength Spectral Plot for a 1D resonator Chain]
    ---------------------------------------------------------
%}

% Objective:
%   Generate a spectral plot for a resonator chain. 
%   A complete derivation of the formula used in this implementation 
%   can be found in arXiv:2409.14537.


% --- Set the parameters ---
    delta = 0.1;        % Contrast (jump of the normal derivative)
    s_1   = 0.6;        % Spacing between the resonators
    N_points = 100;     % Discretisation points used in plot

% --- Define the Bandfunctions in the subwavelength regime ---
    f_1 = @(x) sqrt(delta * (2 / s_1) * (1 - cos(x)));
    f_2 = @(x) sqrt(delta * (2 / s_1) * (1 + cosh(x)));

% --- Define the range for x ---
    x_range = linspace(-pi, pi, N_points);         % First Brillouin zone
    beta_range_2 = linspace(-4, 4, N_points);      % Real numbers

% --- Evaluate the functions on the grid ---
    y_values_1 = real(f_1(x_range));               % Real band 
    y_values_2 = real(f_2(beta_range_2));          % Complex band 

% --- Plotting the spectral Bands ---
    figure;
    lw = 4.5;   % Linewidth
    fs = 18;    % Fontsize
    
    % Plotting of the Bands
    plot(x_range, y_values_1, 'k', 'LineWidth', lw);
    hold on;
    pbaspect([8 6 1]); % Set the ratio to 10:6
    plot(beta_range_2, y_values_2, 'r', 'LineWidth', lw);
    
    % Horizontal lines separating the Band and Bandgap
    y_val = f_1(pi);
    line([-4, 4], [y_val, y_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    
    % Adding labels, legend and ticks
    xlabel('\alpha and \beta respectively', 'FontSize', fs);
    ylabel('\omega^{\alpha, \beta}', 'FontSize', fs);
    legend('Band function', 'Gap function', 'Location', 'northeast', 'FontSize', fs+2);
    set(gca, 'FontSize', fs+4); 
    saveas(gcf, 'Band_1D_Single.pdf');
    hold off;
