%{
    --------------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Find (alpha, beta) such that the SLP has non-trivial kernel]
    ---------------------------------------------------------------------------
%}

clear all;
close all;

% --- Set the parameters ---
    alpha      = [pi, pi];  % Real qusimomentum
    N_multi    = 4;         % Number of multipoles
    N_lattice  = 20;        % Lattice sum size 
    R          = 0.1;       % Resontor radius
    lowerBound = 0;
    upperBound = 5;

% --- compute the minimum in the specified range ---
    options = optimset('TolX', 1e-6, 'TolFun', 1e-6);
    [minValue, minBeta] = fminbnd(@(beta) minSV(beta, alpha, N_multi,N_lattice,  R), lowerBound, upperBound);

% --- Display results ---
    fprintf('----------------------------\n');
    disp(['Minimum SV: ', num2str(minBeta)]);
    disp(['At beta:    ', num2str(minValue)]);

%% --- Double check the result by pltting the smallest singular value ---

% --- plot the smallest singular values ---
x_values = linspace(0, 2*pi, 100); 
plotSmallestSingularValues(x_values, R, alpha, N_multi, N_lattice)


%% --- Definining functions ---


function result = minSV(bet, alp, N_multi, N_lattice, R)

% --- Define the parameters ---
    k0 = 0.001;           
    N = 1;   
    D = 1;          
    c1 = 1/2*D*[1,1];
    c = c1;
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L2 = [L2x,L2y];

    beta  = [bet, bet];
    alpha = alp; 

% --- Generate the single layer potenial ---
    matS  = makeS(k0,R,alpha,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice);
    matR  = makeR(k0, R, alpha, beta, N_multi, N_lattice * 2);
    matSR = (-matR + matS);  

% --- Generat minimal singular values ---
    singularValues = svd(matSR);
    result = min(singularValues);

end


function plotSmallestSingularValues(x_values, R, alpha, N_multi, N_lattice)
   
    % --- Define the parameters ---
        k0 = 0.001;           
        N = 1;   
        D = 1;          
        c1 = 1/2*D*[1,1];
        c = c1;
        d_zeta=makezetadata;
        JHdata = makeJHdata0(k0,R,N_multi);
        JHijdata = makeJHijexpdata(k0,c,N_multi);
        L1x = D;
        L2x = 0;
        L2y = D;
        L2 = [L2x,L2y];

    % --- Loop over the beta range ---
        smallestSingularValues = zeros(size(x_values));
    
        for i = 1:length(x_values)
            x = x_values(i);
            beta = [x, x];
            
            % --- Generate the single layer potenial ---
                matS = makeS(k0, R, alpha, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
                matR = makeR(k0, R, alpha, beta, N_multi, N_lattice * 2);
                matSR = (-matR + matS);
    
            % --- Compute smallest singular value ---
                singularValues = svd(matSR);
                smallestSingularValues(i) = min(singularValues);
        end

    % --- Plot the result ---
        figure;
        plot(x_values, smallestSingularValues, 'o-', 'MarkerFaceColor', 'r');
        xlabel('Value of x');
        ylabel('Smallest Singular Value');
        title('Smallest Singular Values of Matrices Based on Varying x');
        grid on;
end
