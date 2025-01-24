%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [1D skin effect, exponential decay estimate]
    --------------------------------------------------------------
%}

close all;
clear all;

% --- Set the parameters ---
    N_cells     = 20;          % Number of unit cells.
    spacing     = [3, 2, 3];   % Spacings.
    res_length  = [1, 2, 1 ];  % Resonators (same array size as spacings).
    gamma       = 1;           % Complex Gauge potential.

% --- Initialise arrays ---
    s = repmat(spacing,    1, N_cells);  
    l = repmat(res_length, 1, N_cells);   
    n  = length(res_length);          

    N = N_cells * n; % Total number of resonators


%%  --- Plot the eigenvector entries and illustrate exponential decay ---

    plot_all_eigenvectors(N, s, gamma, l, res_length);


%% --- Plot the geometry of the resonator chain ---

    PlotGeometry(N, l, s);


%% --- Plot the eigenvalues ---

    %PlotEigenvalues(N, s, gamma, l)


%% --- Defining functions ---


function plot_all_eigenvectors(N, s, gamma, l, len)

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

        for i = 1:(N)
            plot(1:N, (abs(sortedEigenvectors(:, i))), '.', 'Color', 0.5 * [1, 1, 1], 'LineWidth', lw, 'MarkerSize', 15);
        end

    % --- Add the predicted exponential decay rate ---
        x = 0:N; 
        constant = 1;
        y =  constant * exp(- gamma * 0.5 * sum(len) * (1 / length(len)) * x);

        % Remark:   We are plotting the mode evaluated at the resonators,
        %           therefore the decay rate has to be devided by the number of
        %           resonators in the unit cell, i.e. length(len).

        plot(x, (y), 'r', 'LineWidth', lw); 
   
    % --- Formatting the plot ---
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
                if i == 1
                    capmat(i,j) = (gamma / s(i)) * (ell(i) / (1 - exp(-gamma * ell(i))));
                % Case 2: 1 < i = j < N
                elseif i > 1 && i < N
                    capmat(i,j) = (gamma / s(i))   * (ell(i) / (1 - exp(-gamma * ell(i)))) ...
                                - (gamma / s(i-1)) * (ell(i) / (1 - exp( gamma * ell(i))));
                % Case 3: i = j = N
                else
                    capmat(i,j) = -(gamma / s(i-1)) * (ell(i) / (1 - exp(gamma * ell(i))));
                end

            % --- Populate lower diagonal ---
                elseif i == j - 1
                    % Case 4: 1 <= i = j - 1 <= N - 1
                    capmat(i,j) = -(gamma / s(i)) * (ell(i) / (1 - exp(-gamma * ell(j))));

            % --- Populate upper diagonal ---
                elseif i == j + 1 
                    % Case 5: 2 <= i = j + 1 <= N
                    capmat(i,j) = (gamma / s(j)) * (ell(i) / (1 - exp(gamma * ell(j))));
            end
        end
    end
end

function PlotGeometry(N, l, s)

    x_start = 0;  % Initialize starting position of the first segment            

    figure;
    hold on;
    for i = 1:N
        x_end = x_start + l(i);
        line([x_start, x_end], [0, 0], 'Color', 'k', 'LineWidth', 4);
        if i < N
            x_start = x_end + s(i);
        end
    end

    % --- Formatting the plot ---
    set(gcf, 'Position', [100, 100, 800, 200]);  
    xlabel('x-position', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'FontSize', 18-4);
    set(gca, 'YTick', []); 
    axis equal;  
    grid off;  
    hold off;
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



