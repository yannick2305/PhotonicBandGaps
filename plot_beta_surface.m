% Define the range of beta values
num_points = 100; % Adjust as needed


beta_x = linspace(0, 23, num_points);
beta_y = linspace(0, 23, num_points);
[Beta_X, Beta_Y] = meshgrid(beta_x, beta_y);

W_real = arrayfun(@(x, y) my_function([x, y]), Beta_X, Beta_Y); % Corrected this line


% Create the surface plot
gca = surf(Beta_X, Beta_Y, W_real);

% Remove black lines around tiles
%set(gca, 'EdgeColor', 'none');

grid off;

% Set tick marks on x and y axes to multiples of pi
xticks(0:pi:6*pi);
yticks(0:pi:6*pi);

% Set tick labels to multiples of pi with the symbol "π"
xticklabels({'0','π','2π', '3π', '4π', '5π', '6π'});
yticklabels({'0','π','2π', '3π', '4π', '5π', '6π'});

xlabel('\beta_1', 'FontSize', 22); % Adjust FontSize as needed
ylabel('\beta_2', 'FontSize', 22); % Adjust FontSize as needed
set(gca, 'FontSize', 20); % Adjust FontSize for ticks as needed

% Set figure properties for high resolution
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3]); % Square format
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 3 3]); % Set the size of the figure
set(gcf, 'Renderer', 'painters'); % Choose renderer

% Save as PDF
print('surf_plot_high_2D.pdf', '-dpdf', '-r500'); % 300 dpi resolution


%zlabel('w\_real');


% Define the function f
function w_real = my_function(beta)

    % Set all the Paramters to the same values as in the RUN_band file
    alpha = [pi,pi]; 

    k0 = 0.001;
    N = 1;
    R = 0.1;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = [c1];
    N_lattice = 5;      % Use about 5
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
    CR = makeCR(k0, R, alpha, beta, N_SLP, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));
    w_real = real(ws);   % return only the real part
    upper_bound = 8;
    lower_bound = 0;

    if w_real < lower_bound         % remove the spikes in the singularity
        w_real = lower_bound;
    end 
    if w_real > upper_bound         % remove the spikes in the singularity   
        w_real = upper_bound;
    end

end

