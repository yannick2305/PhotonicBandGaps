clear all;
close all;

%%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
D = 1;              % D = 1 has to be true because of the dual lattice considered
L1x = D;
L2x = 0;
L2y = D;
L1 = [L1x, 0];
L2 = [L2x,L2y];

%%% Radius of bubble
R =  0.05;     %R = 0.05 (Default setting)
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

%%% Set symmetry points in reciprocal space (square lattice)
pointG = [0,0];
pointX = 1/2*b1;
pointM = 1/2*(b1 + b2);

lengthXG = norm(pointG-pointX);
lengthGM = norm(pointM-pointG);
lengthMX = norm(pointX-pointM);
lengthXGM = lengthXG + lengthGM + lengthMX;

%%% Set positions of two bubbles in unit cell
%c1 = 1/2*D*[1,1];
c1 = D*[2/4,1/2];
c2 = D*[3/4,1/2];
N = 2; 
c = [c1;c2];



% Set truncation parameters --------------------
N_multipole = 3;           % use around 7
N_SLP = N_multipole;        % they have to be the same dimension because we are summing them later on        
N_lattice = 3;             % use around 15  
N_XG = 30;                 % 60 (Number of discretization points for beta path)
% --------------------------------------------------------------------

N_GM = floor(sqrt(2)*N_XG);
N_MX = N_XG;

%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;

%%% Plot the geometry
% figure, hold on
% t = linspace(0,2*pi);
% for i = -2:2
%     for j = -2:2
%         for n = 1:N
%             plot(c(n,1)+R*cos(t)+i*L1(1)+j*L2(1), c(n,2)+R*sin(t)+i*L1(2)+j*L2(2),'k')  
%             %text(c(n,1),c(n,2),num2str(n))
%         enda
%     end
% end
% daspect([1 1 1])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set wave numbers ( k0: background, kb: bubble )
v0 = sqrt(kappa_0/rho_0);
vb = sqrt(kappa_b/rho_b);

k0 = 0.001; %0.001
JHdata = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);

% set alpha path 
alps = [(0:N_XG-1)*lengthXG/N_XG, (0:N_GM-1)*lengthGM/N_GM + lengthXG, (0:N_MX-1)*lengthMX/N_MX + lengthXG + lengthGM];

% set beta parth equal to the alpha path
betas = alps;

N_a = length(alps);
bandGM = zeros(N,N_GM);
bandMX = zeros(N,N_MX);
bandXG = zeros(N,N_XG);
    
for alp_ind = 1:N_XG        % Path number 1
    %%% Set quasi-periodic parameter \alpha
    alp =  pointX + (alp_ind-1)/N_XG * (pointG-pointX);
    
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    ws = sort(vb*sqrt(delta*eig(C)./vol));
    bandXG(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_GM     % Path number 2
    %%% Set quasi-periodic parameter \alpha
    alp =  pointG + (alp_ind-1)/N_GM * (pointM-pointG);
    
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    ws = sort(vb*sqrt(delta*eig(C)./vol));
    bandGM(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_MX     % Path number 3
    %%% Set quasi-periodic parameter \alpha
    alp =  pointM + (alp_ind-1)/N_MX * (pointX-pointM);
    
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    ws = sort(vb*sqrt(delta*eig(C)./vol));
    bandMX(:,alp_ind) = ws;
end
band0 = [bandXG, bandGM, bandMX];

%% Plot the complex band __________________________________________________

alphaG = [0, pi];          % fix alpha to the edge of the Brillouin zone.

N_tot = N_GM + N_MX;        % discretization length along triangle.

bandGapXG = zeros(N,N_XG);
bandGapMX = zeros(N,N_tot);

% Beta path parameters----------------------------------
slope  = 0;         % slope of the beta path
length = 2*pi;      % length of beta, i.e. |\beta|
                    % The length above has to amch the tick of the top axis.
% ------------------------------------------------------

length = length / sqrt(1 + slope^2); % scale and normalize length
pointXb = [slope * length, length];

% Essentially we do the two paths instead of just one such that the origin
% of the two plots align, but the scale on the axis may be different.

for alp_ind = 1:N_XG        % Path left to the origin
    %env = 0.5;           
    %pos = pi-(env/2);
    %bet = [pos, pos] + (alp_ind-1)/N_XG *[env,env];
    bet =  -(pointXb) + (alp_ind-1)/N_XG *( pointXb);
    if norm(bet, 1) < 5
        alpha = alphaG;     % fix alpha to the edge of the brillouin zone
   
        C = makeCR_2(k0, R, alpha, bet, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
               
        ws = sort(vb*sqrt(delta*eig(C)./vol));
    else
        ws = NaN;
    end

    bandGapXG(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_tot     % Path right to the origin
    
    bet =  (pointG) + (alp_ind-1)/N_tot * ((1 + sqrt(2))*(pointXb)-pointG);
    
    if norm(bet, 1) < 5
        alpha = alphaG;     % fix alpha (to the edge of the Brillouin zone)
   
        C = makeCR_2(k0, R, alpha, bet, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
        ws = sort(vb*sqrt(delta*eig(C)./vol));
    else
        ws = NaN;
    end

    bandGapMX(:,alp_ind) = ws
   
end

band0G_1 = [bandGapXG, bandGapMX];

bound = 2.48;

bandGapXG_2 = zeros(N,N_XG);
bandGapMX_2 = zeros(N,N_tot);

for alp_ind = 1:N_XG        % Path left to the origin
    %env = 0.5;           
    %pos = pi-(env/2);
    %bet = [pos, pos] + (alp_ind-1)/N_XG *[env,env];
    bet =  -(pointXb) + (alp_ind-1)/N_XG *( pointXb);
    
    if norm(bet, 1) < bound
        alpha = [pi, 0];     % fix alpha to the edge of the brillouin zone
   
        C = makeCR_2(k0, R, alpha, bet, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
               
        ws = sort(vb*sqrt(delta*eig(C)./vol));
    else
        ws = NaN;
    end

    bandGapXG_2(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_tot     % Path right to the origin
    
    bet =  (pointG) + (alp_ind-1)/N_tot * ((1 + sqrt(2))*(pointXb)-pointG);
    
    if norm(bet, 1) < bound
        alpha = [pi, 0];     % fix alpha (to the edge of the Brillouin zone)
   
        C = makeCR_2(k0, R, alpha, bet, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
        ws = sort(vb*sqrt(delta*eig(C)./vol));
    else
        ws = NaN;
    end

    bandGapMX_2(:,alp_ind) = ws
   
end

band0G_2 = [bandGapXG_2, bandGapMX_2];

bandGapXG_3 = zeros(N,N_XG);
bandGapMX_3 = zeros(N,N_tot);

for alp_ind = 1:N_XG        % Path left to the origin
    %env = 0.5;           
    %pos = pi-(env/2);
    %bet = [pos, pos] + (alp_ind-1)/N_XG *[env,env];
    bet =  -(pointXb) + (alp_ind-1)/N_XG *( pointXb);


    if norm(bet, 1) < bound
        alpha = [pi, pi];     % fix alpha to the edge of the brillouin zone
        C = makeCR_2(k0, R, alpha, bet, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);          
        ws = sort(vb*sqrt(delta*eig(C)./vol));
    else
        ws = NaN;
    end

    bandGapXG_3(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_tot     % Path right to the origin
    
    bet =  (pointG) + (alp_ind-1)/N_tot * ((1 + sqrt(2))*(pointXb)-pointG);
    
    if norm(bet, 1) < bound
        alpha = [pi, pi];     % fix alpha (to the edge of the Brillouin zone)
        C = makeCR_2(k0, R, alpha, bet, N_SLP, L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
        ws = sort(vb*sqrt(delta*eig(C)./vol));
    else
        ws = NaN;
    end

    bandGapMX_3(:,alp_ind) = ws
   
end

band0G_3 = [bandGapXG_3, bandGapMX_3];

%% Plots the solution and saves good images ------------------------------
close all
fsz = 16;
lfsz = 18;
alw = 1;
lw = 3;

figure
hold on
w_real   = sort(real(band0),1);
w_real_gap = sort(real(band0G_1),1);
w_imag_gap = abs(sort(imag(band0G_1),1));

w_real_gap_2 = sort(real(band0G_2),1);
w_imag_gap_2 = sort(imag(band0G_2),1);

w_real_gap_3 = sort(real(band0G_3),1);
w_imag_gap_3 = sort(imag(band0G_2),1);

for I_w = 1:N
        plot(alps, w_real(I_w,:), 'b', 'LineWidth', lw)
        %plot(alps, w_real_gap(I_w,:), 'r', 'LineWidth', lw)
        %plot(alps, w_imag_gap(I_w,:), 'g', 'LineWidth', lw)
end
plot(alps, w_real_gap(2,:), 'r', 'LineWidth', lw)
plot(alps, w_real_gap_2(2,:), 'm', 'LineWidth', lw)
plot(alps, w_real_gap_3(1,:), 'c', 'LineWidth', lw)

%plot(alps, w_imag_gap(2,:), 'g', 'LineWidth', lw)
plot(alps, w_imag_gap_2(2,:), 'g', 'LineWidth', lw)
plot(alps, w_imag_gap_3(1,:), 'g', 'LineWidth', lw)


% draw a horizontal line to indicate where the Bandgap starts
hold on
yline(0.505, '--', 'Color', 'k');
yline(0.6955, '--', 'Color', 'k');
hold off

% Tick for the bottom axis
x = [ 0*pi, lengthXG ,  lengthXG+lengthGM,      lengthXGM];
y = {'X', '\Gamma', 'M', 'X'};

% Tick for the top axis
x_b = [ 0, 1/(2 + sqrt(2)), 2*1/(2 + sqrt(2)), 3*1/(2 + sqrt(2))];
y_b = {'-2 \pi', '0' ,'2 \pi', '4\pi'};

axis([0,lengthXGM,0,1.6])
pbaspect([10 6 1]); % Set the ratio to 10:6
legend({'$\omega(\alpha, 0)$', '', '$\omega\left(\left[\begin{array}{c} 0 \\ \pi \end{array}\right], \beta = \frac{t}{\sqrt{1}}\left[\begin{array}{c} 0 \\ 1 \end{array}\right]\right)$', '$\omega\left(\left[\begin{array}{c} \pi \\ 0 \end{array}\right], \beta = \frac{t}{\sqrt{1}}\left[\begin{array}{c} 0 \\ 1 \end{array}\right]\right)$', '$\omega\left(\left[\begin{array}{c} \pi \\ \pi \end{array}\right], \beta = \frac{t}{\sqrt{1}}\left[\begin{array}{c} 0 \\ 1 \end{array}\right]\right)$'}, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 14);

% Set the first set of x-axis ticks and labels at the bottom
set(gca, 'XTick', x);
set(gca, 'XTickLabel', y , 'FontSize', 15);
set(gca, 'XAxisLocation', 'bottom');
set(gca, 'XColor', 'blue');
ylabel(' $\omega (\alpha , \beta)$', 'Interpreter', 'latex' , 'FontSize', 19);
xlabel('$\alpha$ along a path on the Brillouin zone', 'Interpreter', 'latex');

% Set the second set of x-axis ticks and labels at the top
ax1 = gca; % Get current axis
ax2 = axes('Position', ax1.Position- [0, +0.09, 0, 0], 'Color', 'none');
set(ax2, 'XTick', x_b);
set(ax2, 'XTickLabel', y_b);
set(ax2, 'XAxisLocation', 'top');
set(ax2, 'XColor', 'red');
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);

% Hide the y-axis for the second plot
set(ax2, 'YColor', 'none');
set(ax2, 'YTick', [], 'FontSize', 15);

% Hide the top x-axis line
set(ax2, 'XAxis', 'top');

% Set the font size
set(gca, 'FontSize', 16);


% Set the x-axis grid
set(gca, 'XGrid', 'off');

% Add label below the bottom x-axis tick
text(x_b(2), -0.5, 'Custom Label', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');

hold on; % Keep the existing plot while adding the second set of ticks

% Set the positions and labels for the second set of ticks
xticks(x_b);
xticklabels(y_b);

% Set the location and color for the x-axis ticks
set(gca, 'XAxisLocation', 'top');
set(gca, 'XColor', 'red');

%xlabel('Quasi-periodicity $\alpha$ and decay $\beta$','interpreter','latex');
ylabel('Frequency $\omega$','interpreter','latex');
ax = gca;

set(gcf, 'Position', [0, 0, 600, 500])
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex');
l.FontSize = lfsz;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height*0.99];

%saveas(gcf, 'Band_2D_diagonal_beta.pdf');
