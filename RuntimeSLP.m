%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Runtime to generate the single layer potential]
    --------------------------------------------------------------
%}

%% --- Runtime for a single computation ---

clear all;
close all;

% --- Define the parameters ---
    N_multi = 4;   
    k0      = 0.001; 
    alpha   = [3, 3];
    beta    = [0.2, 0.2];

% --- System configuration ---
    N = 1;
    R =  0.03;      
    D = 1;          
    c1 = 1/2*D*[1,1];
    c = c1;
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 0.000001; %1e-3;
    vb = 1;

% --- Compute the runtime for the specified lattice size ---
    Lattice_size = 40; 

% --- Time the runtime ---
    tic;
    % --- Computes the single layer potential ---
    result = makeR(k0, R, alpha, beta, N_multi, Lattice_size);
    % --- UNCOMMNT to compute the capacitance matrix ---
    %result = makeCR(k0, R, alpha, beta, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_inf);
    elapsedTime = toc; 

% --- Display the result ---
    fprintf('---------------------------------\n');
    disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);

%% --- Runtime for a range of different lattice sizes ---

clear all;  
close all;   

% --- Define the parameters ---
    N_multi = 3;   
    k0      = 0.001; 
    alpha   = [3, 3];
    beta    = [0.2, 0.2];

% --- System configuration ---
    N = 1;
    R = 0.03;      
    D = 1;         
    c1 = 1/2*D*[1,1];
    c = [c1];
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x, L2y];
    vol = pi*R^2;
    delta = 0.000001; 
    vb = 1;

% --- Define the range of lattice sizes ---
    N_inf_values = 0:2:200;  

% Initialize an array to store the elapsed times
elapsedTimes = zeros(size(N_inf_values));

% --- Loop through each value of N_inf and measure the runtime ------
for i = 1:length(N_inf_values)

    N_inf = N_inf_values(i);  % Set N_inf to the current value

    clear result;   % Clear previous result
    tic;            % Start timing

    % --- Measure total runtime for the SLP calculation ---
    result = makeR(k0, R, alpha, beta, N_multi, N_inf ); 
    elapsedTimes(i) = toc; 

    % --- Skip the first two iterations ---
    if i == 1 || i == 2
        elapsedTimes(i) = NaN;  
    end

end


%% --- Plot the runtime and fit a quadratic curve ---

% --- Fit a quadratic polynomial (y = a*N^2 + c) ----------------------
    validIdx = ~isnan(elapsedTimes);  % Ignore NaN values for fitting
    N_fit = N_inf_values(validIdx);   % Filter valid N_inf values
    T_fit = elapsedTimes(validIdx);   % Filter valid elapsed times
    
    % Fit only the quadratic term (a*N^2 + c)
    p = polyfit(N_fit.^2, T_fit, 1);  
    a = p(1);  % Coefficient of N^2
    c = p(2);  % Constant term
    
    % Generate the trendline
    trendline = a * N_inf_values.^2 + c;

% --- Plot the runtime as a function of N_inf -----------------------
    lw = 2;    
    fs = 16;    
    
    figure('Position', [100, 100, 1000, 400]); 
    plot(N_inf_values, elapsedTimes, 'x', 'LineWidth', lw, 'MarkerSize', 6);
    hold on;
    plot(N_inf_values, trendline, '--r', 'LineWidth', lw+1);  
    xlabel('Lattice size $N$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Runtime (seconds)', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    grid on;
    legend({'Measured Runtime', 'Quadratic Trendline'}, 'FontSize', fs, 'Interpreter', 'latex');


% --- Display the fitted coefficients (optional) --------------------
    fprintf('----------------------------------------------\n');
    disp('Fitted quadratic trendline (y = a*N^2 + c):');
    disp(['a = ', num2str(a)]);
    disp(['c = ', num2str(c)]);

 print(gcf, 'RuntimeSLP.pdf', '-dpdf', '-bestfit');


%% --- Compare the old runtime to the vectorized version

clear all;  
close all;   

% --- Define the parameters ---
    N_multi = 3;   
    k0      = 0.001; 
    alpha   = [3, 3];
    beta    = [0.2, 0.2];

% --- System configuration ---
    N = 1;
    R = 0.03; %R = 0.05;       
    D = 1;          % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = c1;
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x, L2y];
    vol = pi*R^2;
    delta = 0.000001; %1e-3;
    vb = 1;

    N_inf_values = 3:2:50;  

% --- Initialise arrays ---
    runtime_makeR     = zeros(size(N_inf_values));
    runtime_makeRSlow = zeros(size(N_inf_values));

% --- Loop over the lattice sizes ---
    parfor i = 1:length(N_inf_values)
        N_inf = N_inf_values(i);
        
        % --- Measure runtime for makeR ---
        func_makeR = @() makeR(k0, R, alpha, beta, N_multi, N_inf);
        runtime_makeR(i) = timeit(func_makeR);
        
        % --- Measure runtime for makeRSlow ---
        func_makeRSlow = @() makeRSlow(k0, R, alpha, beta, N_multi, N_inf);
        runtime_makeRSlow(i) = timeit(func_makeRSlow);
    end

%% --- Plot the runtimes --- 
    figure;
    plot(N_inf_values, runtime_makeR, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'makeR');
    hold on;
    plot(N_inf_values, runtime_makeRSlow, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'makeRSlow');
    xlabel('N', 'Interpreter', 'tex');
    ylabel('Runtime (seconds)');
    title('Runtime Comparison of makeR and makeRSlow');
    legend('show');
    grid on;
