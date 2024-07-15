% decay as a function of slope

% Define the
    k0 = 0.001;
    N = 1;
    R = 1;
    D = 1; % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = [c1];
    N_lattice = 15;
    N_SLP = 12;
    N_multipole = N_SLP;
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multipole);
    JHijdata = makeJHijexpdata(k0,c,N_multipole);
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
    alpha =[pi, pi];

% Define the function f(beta)
f = @(beta) sort(vb*sqrt(delta*eig(makeCR(k0, R, alpha, beta, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice))./vol)); % Modify this function as needed

% Define the range of slopes
range = 80;
slope_range = linspace(0, 1, range); 

f_values = zeros(range);

% Compute f(beta) for each slope
for i = 1:range
    slope = slope_range(i);
    beta = 10*(1/sqrt(1 + slope^2))*[slope, 1]; % Define beta vector
    f_values(i) = real(f(beta)); % Compute f(beta)
end

% Plot f(beta) as a function of slope
plot(slope_range, f_values, 'b-', 'LineWidth', 2);
xlabel('Slope');
ylabel('f(\beta)');
title('Plot of f(\beta) as a Function of Slope');
grid on;

