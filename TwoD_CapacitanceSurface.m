%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Singularities of the Capacitance matrix]
    ---------------------------------------------------------
%}

% Objective:
%   For a single resonator inside the unit cell the capacitance matrix
%   is given by a scalar. In this plot, we illustrate the singular behaviour
%   of the capacitance matrix for specific values beta.
%   More details can be found in arXiv:2409.14537.

close all;
clear all;

% --- Set the parameters ---
    alpha = [pi, pi];       % Real quasimomentum
    R     = 0.005;          % Resonator radius
    num_points = 80;        % Resolution of the surface plot

% --- Initialize the surface plot --- 
    % --- Wide view ---
    %beta_x = linspace(0, 20, num_points);
    %beta_y = linspace(0, 20, num_points);

    % --- Close view ---
    beta_x = linspace(2.8, 3.8, num_points);
    beta_y = linspace(2.8, 3.8, num_points);
    
    [Beta_X, Beta_Y] = meshgrid(beta_x, beta_y);

% --- Compute the capacitance matrix at the grid points --- 
    W_real = arrayfun(@(x, y) my_function([x, y], alpha, R), Beta_X, Beta_Y);

% --- Plotting the spectral Bands ---
    figure; 
    ax = axes;  
    fs = 20;    % Fontsize for annotation
    
    % --- Generate the surface plot on the specified axes ---
    surf_handle = surf(ax, Beta_X, Beta_Y, W_real); % Get the surface object handle

    % Remove the black edges around tiles 
    set(surf_handle, 'EdgeColor', 'none'); % Uncomment if necessary

    % Axis labeling
    xticks(ax, 0:pi:6*pi);
    yticks(ax, 0:pi:6*pi);
    zticks(ax, 0:20:20);

    xticklabels(ax, {'$0$', '$\pi$', '$2\pi$', '$3\pi$', '$4\pi$', '$5\pi$', '$6\pi$'});
    yticklabels(ax, {'$0$', '$\pi$', '$2\pi$', '$3\pi$', '$4\pi$', '$5\pi$', '$6\pi$'});
    zticklabels(ax, {'$0$', '$20$'});

    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Set the font size for the tick labels using the Axes handle
    ax.XAxis.FontSize = fs;  
    ax.YAxis.FontSize = fs;  
    ax.ZAxis.FontSize = fs;  

    cb = colorbar;
    cb.FontSize = fs; 
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';

    
    % Add axis labels with custom font size
    xlabel(ax, '$\beta_1$', 'Interpreter', 'latex', 'FontSize', fs);  
    ylabel(ax, '$\beta_2$', 'Interpreter', 'latex', 'FontSize', fs);  
    zlabel(ax, '$\omega^{\alpha, \beta}$', 'Interpreter', 'latex', 'FontSize', fs + 2);
    
    axis equal;
    % --- Set the viewing angle ---
    view(2);         % Top view
    %view(-19, 15);    % Side view
    
    % Save the plot as a PDF
    print('surf_plot_high_2DTop.pdf', '-dpdf', '-r500');

    
%% --- Define the function f ---

function w_real = my_function(beta, alpha, R)

    % Parameters for numerical computations
    N_lattice = 10;      % Use about 10
    N_SLP     = 2;       % Use about 3

    % Parameters/Geometry of the resonator screen
    k0 = 0.001;
    N = 1;
    D = 1;     
    c1 = 1/2*D*[1,1];
    c = c1;
    N_multi = N_SLP;
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L2 = [L2x,L2y];
    vol = pi*R^2;
    kappa_b = 1;
    rho_0 = 1000;
    rho_b = 1;
    delta = rho_b/rho_0;
    vb = sqrt(kappa_b/rho_b);

    % Compute the Capacitance matrix
    CR = makeCR(k0, R, alpha, beta, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    
    % Compute the Bandfunctions
    ws = sort(vb*sqrt(delta*eig(CR)./vol));
    w_real = real(ws);   % return only the real part (imaginary part has to be 0 for admissible solution)
    upper_bound = 8;
    lower_bound = 0;
    
    % Remove the spikes for better visualisation
    if w_real < lower_bound        
        w_real = lower_bound;
    end 
    if w_real > upper_bound         
        w_real = upper_bound;
    end

end
