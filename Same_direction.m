
% Suppose that alpha and beta have the same direction
% Define the slope of this direction


num_points = 100; % Adjust as needed

alpha = linspace(0, 2*pi +0.1, num_points);
beta  = linspace(0, pi, num_points);
[Alpha, Beta] = meshgrid(alpha, beta);

W_real = arrayfun(@(x, y) my_function([x, y]), Alpha, Beta); 

% Create the surface plot
gca = surf(Alpha, Beta, W_real);

% Set tick marks on x and y axes to multiples of pi
xticks(0:pi:5*pi);
yticks(0:pi:5*pi);

% Set tick labels to multiples of pi with the symbol "π"
xticklabels({'0','π','2π', '3π', '4π', '5π'});
yticklabels({'0','π','2π', '3π', '4π', '5π'});

xlabel('\alpha', 'FontSize', 22); % Adjust FontSize as needed
ylabel('\beta', 'FontSize', 22); % Adjust FontSize as needed
%set(gca, 'FontSize', 20); % Adjust FontSize for ticks as needed

% Define the function f
function w_real = my_function(alphabeta)

    alpha = alphabeta(1);
    beta = alphabeta(2);

    slope = 1;

    alp = [alpha * slope, alpha]; %alpha and beta along the same slope
    bet = [beta * slope, beta];

    % Set all the Paramters to the same values as in the RUN_band file
    k0 = 0.001;
    N = 1;
    R = 0.05;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = [c1];
    N_lattice = 4;      % Use about 5
    N_SLP = 4;          % Use about 5
    N_multi = N_SLP;
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    kappa_b = 1;
    kappa_0 = 1000;
    rho_0 = 1000;
    rho_b = 1;
    delta = rho_b/rho_0;
    v0 = sqrt(kappa_0/rho_0);
    vb = sqrt(kappa_b/rho_b);

    % Compute the Bandfunctions
    CR = makeCR(k0, R, alp, bet, N_SLP, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));
    w_real = real(ws);   % return only the real part
    w_imag = imag(ws);

    tol = 2;
    if w_imag > tol        % remove the spikes in the singularity
        w_real = 0;
    end 

    upper_bound = 5;
    lower_bound = 0;
 
    if w_real < lower_bound        % remove the spikes in the singularity
        w_real = lower_bound;
    end 
    if w_real > upper_bound       % remove the spikes in the singularity   
        w_real = upper_bound;
    end



end
