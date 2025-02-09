%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Sorface plot for nontrivial Kernel of the SLP]
    --------------------------------------------------------------
%}

%%
clear all;
close all;

% --- Fix the parameters ---
    N_multi    = 4;         % Use about 4 
    N_lattice  = 5;         % Use at least 5
    alpha      = [pi, pi]; 
    Resolution = 200;

% --- Define the parameters ---
    k0 = 0.0001;
    N = 1;
    R = 0.05;       
    D = 1;         
    c1 = 1/2*D*[1,1];
    c = c1;
    d_zeta = makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 0.001; 
    vb = 1;

% --- Precompute the real SLP once ---
    matS = makeS(k0,R,alpha,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice);
    singularValues_S = svd(matS);
    minsing = min(singularValues_S);

% --- Large View ---
    %beta1_values = linspace(0, 4 * pi, Resolution); 
    %beta2_values = linspace(0, 4 * pi, Resolution); 

% --- Close up on the first nontrivil kernel ---
    beta1_values = linspace(2.9, 4.05, Resolution); 
    beta2_values = linspace(2.9, 4.05, Resolution); 

% --- Initialize an array ---
    smallestSingularValues = zeros(length(beta1_values), length(beta2_values));

% --- Loop over the beta grid ---
    for i = 1:length(beta1_values)
        for j = 1:length(beta2_values)
            beta1 = beta1_values(i);
            beta2 = beta2_values(j);
            beta = [beta1, beta2];
            
            % Calculate matR matrices
            matR = makeR(k0, R, alpha, beta, N_multi, N_lattice);
            
            % Calculate matSR and its smallest singular value
            matSR = - matR + matS; 
            singularValues = svd(matSR);
            smallestSingularValues(i, j) = min(min(singularValues),0.1);
        end
    end
    
    minValue = min(smallestSingularValues(:));

% --- Plot the solution ---
    [Beta1, Beta2] = meshgrid(beta1_values, beta2_values);
    
    figure;
    surf(Beta1, Beta2, smallestSingularValues, 'EdgeColor', 'none');
    hold on;
    
    % Add labels and title
    xlabel('$\beta_1$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\beta_2$', 'FontSize', 18, 'Interpreter', 'latex');
    xticks(0:pi:3*pi);  
    yticks(0:pi:3*pi);   
    
    % Customize x-axis labels to display multiples of pi
    xticklabels({'$0$', '$\pi$', '$2\pi$', '$3\pi$'});
    yticklabels({'$0$', '$\pi$', '$2\pi$', '$3\pi$'});
    set(gca, 'TickLabelInterpreter', 'latex');
    
    % Increase font size of tick labels
    ax = gca; 
    ax.XAxis.FontSize = 22; 
    ax.YAxis.FontSize = 22; 
    ax.ZAxis.FontSize = 22; 

    % Add a color bar and enable grid
    colorbar;
    axis tight; 
    axis equal;
    grid on;
    
    % Get the current colorbar object
    cb = colorbar;
    
    % Change the font size of the colorbar labels
    cb.FontSize = 18; 
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    view(0, 90);
    
    hold off;
    
