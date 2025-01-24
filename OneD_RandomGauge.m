%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [1D skin effect with a random gauge potential]
    --------------------------------------------------------------
%}

close all;
clear all;

% --- Set the parameters ---
    N_cells = 40;          % Number of unit cells.
    n       = 2;           % Number of resonators per unit cell.
    spacing = [1, 2];      % Spacings   (length n).
    length  = [0.8, 1.2];  % Resonators (length n).

    gamma_mean = 1;        % Mean of the gauge potential.
    gamma_min  = 0.5;      % Lower bound of gauge potential.
    gamma_max  = 2;        % Upper bound of gauge potential.


% --- Generate random gauge potential on resonators ---

    M = N_cells * n;       % Length of the resonator chain
    
    % --- Generate random values in the range ---
    array = gamma_min + (gamma_max - gamma_min) * rand(1, M); 
    
    % --- Adjust to the desired mean ---
    current_mean = mean(array);
    array = array - current_mean + gamma_mean; 

    % --- Clip values to ensure they stay in the range ---
    array = max(gamma_min, min(gamma_max, array));
    
    % --- Re-normalize to keep the mean ---
    final_mean = mean(array);
    array = array + (gamma_mean - final_mean); 
    array = max(gamma_min, min(gamma_max, array)); 

    % --- Display the generated array ---
    fprintf('----------------------------------\n');
    fprintf('Mean: %.2f, Min: %.2f, Max: %.2f\n', mean(array), min(array), max(array));

    gamma = array;

% --- Initialise resonator chain ---
    s = repmat(spacing, 1, N_cells);    % Periodic spacings.
    l = repmat(length , 1, N_cells);    % Periodic lengths.
    N = N_cells * n;                    % Total number of resonators.


%%  --- Plot the eigenvector entries and illustrate exponential decay ---

    plot_all_eigenvectors(N, s, gamma, l, length, gamma_mean);


%% --- Plot the geometry of the resonator chain ---

    PlotGeometry(N, l, s, gamma);


%% --- Plot the eigenvalues ---

    % PlotEigenvalues(N, s, gamma, l)


%% --- Defining functions ---

function plot_all_eigenvectors(N, s, gamma, l, len, mean)

    lw = 3;     % Thickness of the resonators
    fs = 24;    % Fontsize of the annotations

    % --- Generate the eigenvectors of the Capacitance matrix ---
        capmat = Capacitance(N, s, gamma, l);
        [V, D] = eig(capmat);
        eigenvalues = diag(D);
        [~, sortIdx] = sort(eigenvalues);
        sortedEigenvectors = V(:, sortIdx);    
    
    % --- Plot the eigenvector entries ---
        figure;
        hold on;    

        for i = 2:(N)
            plot(1:N, (abs(sortedEigenvectors(:, i))), '-', 'Color', 0.5 * [1, 1, 1], 'LineWidth', lw, 'MarkerSize', 15);
        end

    % --- Add the predicted exponential decay rate ---
        x = 0:N; 
        constant = 1;
        y =  constant * exp(- 0.5 * (mean) * sum(len) / length(len) * x); 

        % Remark:   We are plotting the mode evaluated at the resonators,
        %           therefore the decay rate has to be devided by the number of
        %           resonators in the unit cell, i.e. length(len).

        plot(x, (y), 'r', 'LineWidth', lw); 
   
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


function capmat = Capacitance(N, s, gamma,  ell)
    
    capmat = zeros(N, N);  
    
    % --- Popolate the matrix ---
    for i = 1:N
        for j = 1:N
            if i == j
            % --- Populate  diagonal ---
                % Case 1: 1 = i = j
                if i == 1
                    capmat(i,j) = (gamma(i) / s(i)) * (ell(i) / (1 - exp(-gamma(i) * ell(i))));
                % Case 2: 1 < i = j < N
                elseif i > 1 && i < N
                    capmat(i,j) = (gamma(i) / s(i))   * (ell(i) / (1 - exp(-gamma(i) * ell(i)))) ...
                                - (gamma(i) / s(i-1)) * (ell(i) / (1 - exp( gamma(i) * ell(i))));
                % Case 3: i = j = N
                else
                    capmat(i,j) = -(gamma(i) / s(i-1)) * (ell(i) / (1 - exp(gamma(i) * ell(i))));
                end

            % --- Populate lower diagonal ---
                elseif i == j - 1
                    % Case 4: 1 <= i = j - 1 <= N - 1
                    capmat(i,j) = -(gamma(j) / s(i)) * (ell(i) / (1 - exp(-gamma(j) * ell(j))));

            % --- Populate upper diagonal ---
                elseif i == j + 1 
                    % Case 5: 2 <= i = j + 1 <= N
                    capmat(i,j) = (gamma(j) / s(j)) * (ell(i) / (1 - exp(gamma(j) * ell(j))));
            end
        end
    end
end


function PlotGeometry(N, l, s, gamma)

    % --- Normalize gamma for color scaling ---
    gamma_min = min(gamma);
    gamma_max = max(gamma);
    normalized_gamma = (gamma - gamma_min) / (gamma_max - gamma_min);   % Normalize to [0, 1]
    cmap = [linspace(0, 1, 256)', zeros(256, 1), linspace(1, 0, 256)']; % Gradient from blue to red

    x_start = 0;     % Starting position         

    figure;
    hold on;

    % --- Draw and color each resonator ---
    for i = 1:N
        % Determine color based on gamma value
        color_index = round(normalized_gamma(i) * 255) + 1; % Scale to [1, 256]
        color_index = max(1, min(256, color_index)); 
        color = cmap(color_index, :);

        % Compute the end position of the resonator
        x_end = x_start + l(i);

        % Plot the resonator with the assigned color
        line([x_start, x_end], [0, 0], 'Color', color, 'LineWidth', 4);

        % Update x_start for the next segment (include spacing)
        if i < N
            x_start = x_end + s(i);
        end
    end

    % --- Formatting the plot ---
    set(gcf, 'Position', [100, 100, 600, 120]);  
    xlabel('x-position', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'FontSize', 14);
    set(gca, 'YTick', []); 
    axis equal;  
    grid off;  
    hold off;

    %  ---Colorbar settings ---
    colormap(cmap);
    c = colorbar;
    c.Ticks = [0, 0.5, 1];
    c.TickLabels = {sprintf('%.2f', gamma_min), '0', sprintf('%.2f', gamma_max)};
    c.Label.String = '$\gamma$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 18;
    saveas(gcf, 'RandomGaugeChainN.pdf');
end


function PlotEigenvalues(N, s, gamma, l)
    
    % --- Generate the Capacitance matrix ---
    capmat = Capacitance(N, s, gamma, l);
    
    % --- Compute the eigenvalues of the capacitance matrix ---
    [~, eigenvalues_matrix] = eig(capmat);
    eigenvalues = sort(diag(eigenvalues_matrix));  
    
    % --- Plot the eigenvalues ---
    figure;
    plot(eigenvalues, 'bo', 'LineWidth', 2, 'MarkerSize', 8);  
    xlabel('Index', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', 14);
    grid on; 
end
