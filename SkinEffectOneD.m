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
    N      = 62;                    % Number of resonators (even number for dimers).
    gamma  = 1;                     % Complex Gauge potential.
    s      = mod(1:(N-1), 2) + 1;   % Resonator spacings (1, 2, 1, 2,..., 1, 2).
    l      = mod(1:N    , 1) + 1;   % Length of the resonators.


%%  --- Plot the eigenvector entries and illustrate exponential decay ---

    plot_all_eigenvectors(N, s, gamma, l);


%% --- Plot the geometry of the resonator chain ---

% --- Plot the resonators ---
    % Initialize positions
    x_start = 0;              
    y_positions = 0;           

    figure;
    hold on;
    for i = 1:N
        x_end = x_start + l(i);
        line([x_start, x_end], [y_positions, y_positions], 'Color', 'k', 'LineWidth', 4);
    
        % Update the starting position for the next resonator
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


%% --- Plot the eigenvalues ---

% --- Generate the Capacitance matrix
    capmat = generate_capmat(N, s, gamma, l);

% --- Compute the eigenvalues of the capacitance matrix ---
    [eigenvectors, eigenvalues_matrix] = eig(capmat);
    eigenvalues = sort(diag(eigenvalues_matrix));

% --- Plot the eigenvalues ---
    figure;
    plot(eigenvalues, 'bo', 'LineWidth', 2, 'MarkerSize', 8);  
    xlabel('Index',      'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;


%% --- Defining functions ---


function plot_all_eigenvectors(N, s, gamma, l)

    lw = 3;     % Thickness of the resonators
    fs = 24;    % Fontsize of the annotations

    % --- Generate the eigenvectors of the Capacitance matrix ---
        capmat = generate_capmat(N, s, gamma, l);
        [V, D] = eig(capmat);
        eigenvalues = diag(D);
        [~, sortIdx] = sort(eigenvalues);
        sortedEigenvectors = V(:, sortIdx); % Reorder the eigenvectors

        figure;
        hold on;    
    
    % --- Plot the eigenvector entries ---
        figure;
        hold on;    

        for i = 1:N
            plot(1:N, (abs(sortedEigenvectors(:, i))), '.', 'Color', 0.5 * [1, 1, 1], 'LineWidth', lw, 'MarkerSize', 15);
        end

    % --- Add the predicted exponential decay rate ---
        x = 0:N; 
        y = exp(- gamma * 0.5 * 0.5 * (l(1)+l(2)) * x); 
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



function capmat = generate_capmat(N, s, gamma, l)
    
    capmat = zeros(N, N);

    % --- Populate the Capacitance matrix ---
    for i = 1:N
        for j = 1:N
            if i == 1 && j == 1
                % Case 1: i = j = 1
                capmat(i,j) = (gamma / s(i)) * (l(i) / (1 - exp(-gamma * l(i))));
            elseif i == j && i > 1 && i < N
                % Case 2: 1 < i = j < N
                capmat(i,j) = (gamma / s(i)) * (l(i) / (1 - exp(-gamma * l(i)))) ...
                    - (gamma / s(i-1)) * (l(i) / (1 - exp(gamma * l(i))));
            elseif i == j - 1
                % Case 3: i = j - 1
                capmat(i,j) = - (gamma / s(i)) * (l(j) / (1 - exp(-gamma * l(j))));
            elseif i == j + 1
                % Case 4: i = j + 1
                capmat(i,j) = (gamma / s(j)) * (l(i) / (1 - exp(gamma * l(j))));
            elseif i == N && j == N
                % Case 5: i = j = N
                capmat(i,j) = - (gamma / s(N-1)) * (l(N) / (1 - exp(gamma * l(N))));
            end
        end
    end
end
