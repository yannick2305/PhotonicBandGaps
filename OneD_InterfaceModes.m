%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [January 2025]
    Description:  [1D Exponentially localised interface modes]
    --------------------------------------------------------------
%}

% Objective: Code supporting Section 2.2.5. of the paper arXiv:2409.14537


clear all;
close all;

% --- Set system parameters ---
    n = 101;
    s_1 = 1;                % Spacing between first and second resonator
    s_2 = 2;                % Spacing between second and first resonator
    delta = 0.001;          % Contrast
    N_points = 200;         % Discretisation points used in plot
    L =  s_1 + s_2 + 2;     % Length of the unit cell.

% --- Get resonant modes and frequencies via the capacitance matrix ---
    C  = generate_matrix(n, s_1, s_2);
    ws = sort(1 * sqrt(delta*eig(C)));

% --- Define the Bandfunctions in the subwavelength regime ---
    f_1 = @(x) sqrt(delta) * sqrt((1/s_1 + 1/s_2)     + sqrt(1/(s_1)^2     + 1/(s_2)^2 + 2/(s_1*s_2)* cos(x * L)));
    f_2 = @(x) sqrt(delta) * sqrt((1/s_1 + 1/s_2)     - sqrt(1/(s_1)^2     + 1/(s_2)^2 + 2/(s_1*s_2)* cos(x * L)));
    f_3 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) - sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 - 2/(s_1*s_2)*cosh(x * L)))));
    f_4 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) + sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 - 2/(s_1*s_2)*cosh(x * L)))));
    f_5 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) + sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 + 2/(s_1*s_2)*cosh(x * L)))));

% --- Define the range for x ---
    x_range = linspace(-pi/L, pi/L, N_points);  % First Brillouin zone.    

% --- Plot the complex bands over some ranges ---
    beta_bound = acosh( (1/2) * (s_2^2 + s_1^2)/(s_1*s_2) ); 
    beta_range = linspace(- (1/L) * beta_bound, (1/L) * beta_bound, N_points);
    beta_range_2 = linspace(-4/L - 0.2, 4/L + 0.2, N_points);   

% --- Evaluate the functions on the grid ---
    y_values_1 = real(f_1(x_range));        % Real Band
    y_values_2 = real(f_2(x_range));        % Real Band
    y_values_3 = real(f_3(beta_range));     % Complex Band
    y_values_4 = real(f_4(beta_range));     % Complex Band
    y_values_5 = real(f_5(beta_range_2));   % Complex Band


% --- Overlay resonant frequencies onto th bandplot ---
    x_search = linspace(-pi/L , pi/L , 1000);
    intersection_points = NaN(length(ws), 1);
    tol = 1e-3;
    
    for i = 1:floor(length(ws)/2+1)
        w = ws(i);
    
        diff_f1 = abs(f_2(x_search) - w);
        [min_diff, min_index] = min(diff_f1);
        
        if min_diff < tol
            intersection_points(i) = x_search(min_index);
        end
    end
    
    for i = floor(length(ws)/2 + 2) : n
        w = ws(i);
    
        diff_f1 = abs(f_1(x_search) - w);
        [min_diff, min_index] = min(diff_f1);
    
        if min_diff < tol
            intersection_points(i) = x_search(min_index);
        end
    end

    w = ws(floor(n/2)+1);
    diff_f1 = abs(f_3(x_search) - w);
    [min_diff, min_index] = min(diff_f1);
    
    gap_decay = x_search(min_index);


% --- Plotting the spectral Bands and resonant frequencies ---
    figure;
    lw = 4.5;   % Linewidth
    fs = 14;    % Fontsize
    
    % Plotting of the Bands
    plot(x_range, y_values_1, 'k', 'LineWidth', lw);

    hold on;
    pbaspect([10 4 1]);  
    
    % Plot the lines
    plot(x_range, y_values_2, 'k', 'LineWidth', lw);
    plot(beta_range, y_values_3, 'r', 'LineWidth', lw);
    plot(beta_range, y_values_4, 'r', 'LineWidth', lw);
    plot(beta_range_2, y_values_5, 'r', 'LineWidth', lw);
    
    % Horizontal lines separating the Band and Bandgap
    y_val = f_1(pi/L);
    line([-4/L, 4/L], [y_val, y_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    b_val = f_2(pi/L);
    line([-4/L, 4/L], [b_val, b_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    g_val = f_1(0);
    line([-4/L, 4/L], [g_val, g_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
    
    % Plot frequencies inside the bandplot
    band_resonance_handles = [];
    for i = 1:length(intersection_points)
        if ~isnan(intersection_points(i))
            h1 = plot(intersection_points(i), ws(i), '|', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 1.2, 'DisplayName', 'Band resonance');
            band_resonance_handles = [band_resonance_handles, h1];  
            h2 = plot(-intersection_points(i), ws(i), '|', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 1.2, 'DisplayName', 'Band resonance');
            band_resonance_handles = [band_resonance_handles, h2];  
        end
    end
  
    gap_resonance_handles = [];
    h3 = plot(-gap_decay, ws(floor(n/2)+1), 'x', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 4, 'DisplayName', 'Gap resonance');
    gap_resonance_handles = [gap_resonance_handles, h3];  
    h4 = plot(gap_decay, ws(floor(n/2)+1), 'x', 'MarkerEdgeColor', 'b', 'MarkerSize', 10, 'LineWidth', 4, 'DisplayName', 'Gap resonance');
    gap_resonance_handles = [gap_resonance_handles, h4];  
    
    plot(0, 0, '|', 'MarkerEdgeColor', 'b', 'MarkerSize', 15,'LineWidth', 4, 'DisplayName', sprintf('w = %.2f', ws(i))); 
    hold off;
    
    % Figure properties
    set(gcf, 'Position', [0, 0, 600, 500]);
    xticks([-pi/L, 0, pi/L]);
    xticklabels({'$-\pi/L$', '$0$', '$\pi/L$'});
    xlabel('$\alpha$ and $\beta$ respectively', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('\omega^{\alpha, \beta}', 'FontSize', fs);
    ylim([0, 0.09]);
    xlim([-pi/L-0.1, pi/L + 0.1]);
    legend([band_resonance_handles(1), gap_resonance_handles(1)], {'Band Resonance', 'Gap Resonance'}, 'Location', 'northeast', 'FontSize', fs+2);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');

    saveas(gcf, 'InterfaceGapFlatNew.pdf');
     
%% --- Plot the eigenmodes ---

    % --- Compute the decay length ---
    defect_freq  = compute_omega_i(delta, s_1, s_2);
    decay_length = compute_beta(L, s_1, s_2, defect_freq, delta);
    
    % --- Plot the eigenmodes and the predicted decay length ---
    shift = 2 ;  % shift the decay envelope if needed
    plot_all_eigenvectors(C, n, decay_length, shift);


%% --- Defining functions ---


function C = generate_matrix(n, s1, s2)
    % Objective: compute the interface capacitance matrix

    if n < 3
        error('Matrix size n must be at least 3.');
    end

    % Define the parameters based on the input equations
    beta1 = -1 / s1;          
    beta2 = -1 / s2;         
    alpha = 1 / s1 + 1 / s2;  
    eta = 2 / s2;             
    tilde_alpha = 1 / s1;     

    C = zeros(n, n);

    % Populate diagonal
    for i = 1:n
        C(i, i) = alpha;
    end

    % Populate upper and lower diagonal
    for i = 1:n-1
        if mod(i, 2) == 1
            C(i, i+1) = beta1; % Upper diagonal
            C(i+1, i) = beta1; % Lower diagonal
        else
            C(i, i+1) = beta2; % Upper diagonal
            C(i+1, i) = beta2; % Lower diagonal
        end
    end

    % Manually set the corners
    C(1, 1) = tilde_alpha;
    C(n, n) = tilde_alpha;

    % Set the defect in the middle
    if mod(n, 2) == 1
        mid = ceil(n / 2);
        C(mid, mid) = eta;
    else
        mid1 = n / 2;
        mid2 = mid1 + 1;
        C(mid1, mid1) = eta;
        C(mid2, mid2) = eta;
    end
end

function plot_all_eigenvectors(C, N, beta, shift)
    % Objective: Plot the eigenvectors of the Capacitance matrix and
    % overlay the predicted exponential decay rate.

    lw = 3;     
    fs = 24;    

    % --- Generate the eigenvectors of the Capacitance matrix ---
        capmat = C;
        [V, D] = eig(capmat);
        eigenvalues = diag(D);
        [~, sortIdx] = sort(eigenvalues);
        sortedEigenvectors = V(:, sortIdx);    
    
    % --- Plot the eigenvector entries ---
        figure;
        hold on;    
        Pos = floor(N/2) + 2;
        for i = Pos  : Pos
            plot(1:N, (abs(sortedEigenvectors(:, i))), '-', 'Color', 0.5 * [1, 1, 1], 'LineWidth', lw, 'MarkerSize', 15);
        end

    % --- Add the predicted exponential decay rate ---
        x = floor(N/2)+shift : N; 
        x2 = 1 : floor(N/2)+shift; 
        constant = exp(-0.8);
        y =  constant * exp(- beta * (x - floor(N/2) - shift)); 
        plot(x, (y), 'r', 'LineWidth', lw); 

        y1 =  constant * exp( beta * (x2 - floor(N/2) - shift)); 
        plot(x2, (y1), 'r', 'LineWidth', lw);
   
    % --- Formatting the plot ---
        set(gcf, 'Position', [100, 100, 400, 400]); 
        xlabel('Eigenvector index', 'Interpreter', 'latex', 'FontSize', fs);
        ylabel('$|u(x)|$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([-25, 1]);
        xlim([0, floor(N/10)*10]);
        grid on;
        set(gca, 'YScale', 'log');
        set(gca, 'FontSize', fs-4);
        hold off;
end

function omega_i = compute_omega_i(delta, s1, s2)
    % Objective: Compute defect frequency
    
    sqrt_terms   = sqrt(9 / s1^2 - 14 / (s1 * s2) + 9 / s2^2);
    linear_terms = (3 / s1) + (3 / s2);
    
    omega_i = sqrt(delta * 0.5 * (-sqrt_terms + linear_terms));
end

function beta = compute_beta(L, s1, s2, omega_i, delta)
    % Objective: Compute the decay length based on the frequency
   
    inner_arcosh = (s1 * s2 / 2) * (1/s1^2 + 1/s2^2 - (1/s1 + 1/s2 - (omega_i^2 / delta))^2);
    beta = (1 / L) * acosh(inner_arcosh);
end
