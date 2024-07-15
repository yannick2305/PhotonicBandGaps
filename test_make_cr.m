%% Test the function parameters: (debugging purposes)

 k = 0.001;             % Set a resonant frequency
 R = 1;                % Set the radius of the resonator
 alpha = [pi, pi];     % Fix alpha to the edge of the Brillouin zone
 beta = [3*pi, 3*pi];  % Beta will be varied later on in the plot
 num_points = 20;      % Discretization in the quadrature
 N_SLP = 5;            % The size of the truncated Fourier Basis
 N_multi = N_SLP;
 N_lattice = 7;       % Size of the truncated lattice
 D = 1; % D = 1 has to be true
 N = 1;
L1x = D;
L2x = 0;
L2y = D;
L1 = [L1x, 0];
L2 = [L2x,L2y];

%%% Radius of bubble
R =  0.05;  %R = 0.05
vol = pi*R^2;

%%% Material parameters for bubble
rho_0 = 1000;
rho_b = 1;
kappa_0 = 1000;
kappa_b = 1;

%%% High contrast parameter \delta
delta = rho_b/rho_0;

%%% Set reciprocal vectors
b1 = 2*pi/(L1x*L2y)*[L2y,-L2x];
b2 = 2*pi/L2y*[0,1];
JHdata = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);

d_zeta=makezetadata;

 %Call the makeR function with the provided variables
 SLP_result = makeCR_l(k, R, alpha, beta, num_points, N_SLP, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);

 %Display or use the result as needed
 disp(SLP_result);

% Define the range of t values
tolerance = 0.05;
t_range = linspace(3*pi-tolerance, 3*pi , 200); % Adjust the range and number of points as needed

% Initialize an array to store the function outputs
output_values = zeros(size(t_range));

% Loop through each t value, create a vector input for the function, evaluate the function, and store the result
for i = 1: 100
    t = t_range(i);
    beta = [t, t]; % Create a vector from the scalar t
    output_values(i) = makeCR_l(k, R, alpha, beta, num_points, N_SLP, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
end

% Plot the results
plot(t_range, real(output_values));
xlabel('t');
ylabel('Function Output');
title('Plot of your function along the line beta = t(1,1)');
