%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Subwavelength Spectral Plot for a 1D Dimer Chain]
    ---------------------------------------------------------
%}

% Objective:
%   Generate a spectral plot for a dimer resonator chain. 
%   A complete derivation of the formula used in this implementation 
%   can be found in arXiv:2409.14537.

clear all;
close all;

% --- Set the parameters ---
    delta = 0.1;        % Contrast (jump of the normal derivative)
    s_1 = 0.8;          % Spacing between first and second resonator
    s_2 = 2;            % Spacing between second and first resonator
    N_points = 100;     % Discretisation points used in plot

% --- Define the Bandfunctions in the subwavelength regime ---
    f_1 = @(x) sqrt(delta) * sqrt((1/s_1 + 1/s_2)     + sqrt(1/(s_1)^2     + 1/(s_2)^2 + 2/(s_1*s_2)*cos(x)));
    f_2 = @(x) sqrt(delta) * sqrt((1/s_1 + 1/s_2)     - sqrt(1/(s_1)^2     + 1/(s_2)^2 + 2/(s_1*s_2)*cos(x)));
    f_3 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) - sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 - 2/(s_1*s_2)*cosh(x)))));
    f_4 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) + sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 - 2/(s_1*s_2)*cosh(x)))));
    f_5 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) + sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 + 2/(s_1*s_2)*cosh(x)))));

% --- Define the range for x ---
    % Plot real Bands over the first Brillouin zone
    x_range = linspace(-pi, pi, N_points);      

    % Plot the complex bands over some ranges
    beta_bound = acosh((1/2) * (s_2^2 + s_1^2)/(s_1*s_2)); 
    beta_range = linspace(-beta_bound, beta_bound, N_points);
    beta_range_2 = linspace(-4, 4, N_points);   

% --- Evaluate the functions on the grid ---
    y_values_1 = real(f_1(x_range));        % Real Band
    y_values_2 = real(f_2(x_range));        % Real Band
    y_values_3 = real(f_3(beta_range));     % Complex Band
    y_values_4 = real(f_4(beta_range));     % Complex Band
    y_values_5 = real(f_5(beta_range_2));   % Complex Band


% --- Plotting the spectral Bands ---
    figure;
    lw = 4.5;   % Linewidth
    fs = 18;    % Fontsize
    
    % Plotting of the Bands
    plot(x_range, y_values_1, 'k', 'LineWidth', lw);
    hold on;
    pbaspect([8 6 1]); % Set the ratio to 10:6
    plot(x_range, y_values_2     , 'k', 'LineWidth', lw);
    plot(beta_range, y_values_3  , 'r', 'LineWidth', lw);
    plot(beta_range, y_values_4  , 'r', 'LineWidth', lw);
    plot(beta_range_2, y_values_5, 'r', 'LineWidth', lw);
    
    % Horizontal lines separating the Band and Bandgap
    y_val = f_1(pi);
    line([-4, 4], [y_val, y_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    b_val = f_2(pi);
    line([-4, 4], [b_val, b_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    g_val = f_1(0);
    line([-4, 4], [g_val, g_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    
    % Adding labels, legend and ticks
    xlabel('\alpha and \beta respectively', 'FontSize', fs);
    ylabel('\omega^{\alpha, \beta}', 'FontSize', fs);
    legend('Band function', '', 'Gap function', '', '', 'Location', 'northeast', 'FontSize', fs+2);
    set(gca, 'FontSize', fs+4); 
    saveas(gcf, 'Band_1D_Dimer.pdf');
    hold off;

