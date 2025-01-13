%{
    --------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Three dimensional resonator chains.]
    ---------------------------------------------------------------------------
%}

% Acknowledgement:
%           This File is just an application of the three dimensional
%           capacitance matrix, originally implemented by Thea KOSHE.

clear all;
close all;

% --- System Geometry ---
    latVect = [1; 0; 0];    % Arrangement of the unit cells.
    numfund = 400;          % Number of unit cells.
    posRes  = [0.5; 0; 0];  % Resonator position inside unit cell.
    numRes  = 1;            % Number of resonators in the unit cell.
    radRes  = 0.01;       % Resonator radius.
    N_mul   = 0;            % Multipole expansion.
    
% --- Define the metamaterial ---
    infMat = infMetaMat(latVect,numRes,posRes,radRes);

% --- Plot the metamaterial ---
    infMat.plot(4)

% --- Compute the capacitance matrix ---
    [~,finMat] = finSingleLayerMatrix(infMat,N_mul,numfund);
    [~,capMat] = capMatrix(finMat);

% --- Computing frequencies (to find Bandgap) ---
    [V,D]= eig(capMat);
    frequencies = sqrt(diag(D));


%% --- Compute the Decay ---

% --- Initialise ---
    C   = capMat; 
    I   = eye(numfund * numRes);
    Nx  = numfund * numRes;
    w   = max(frequencies) + 0.5;   % Chose a frequency inside the bandgap.

% --- Compute the discrete Greens function ---
    E  = C - w^2 * I;
    In = E \ I;     

    i =  floor(numfund / 2);    % Avoid edge effects.
    column_i = In(:, i);        

    row = column_i(floor(numfund / 2) : end);  

% --- Plot the vector entries ---
    figure;

    loglog(1:length(row), abs(row), '-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8);

    xlabel('$|i-j|$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$C_{ij}$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'FontSize', 14);
    pbaspect([10 3 1]);
    grid on;

% --- Compute the algebraic decay rate ---

    % Remove possible edge effects
    disf = 2 ;
    disl = floor(numfund/4);
    x = (1 +disf):(length(row) - disl);  
    y = abs(row((1 +disf):end-disl)); 

    % Fit a linear curve
    log_x = log(x);
    log_y = log(y);
    coefficients = polyfit(log_x, log_y, 1);
    slope = coefficients(1);
    
    fprintf('-------------------------------------\n');
    fprintf('The decay rate is: O(n^{%.4f})\n', slope);
