%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [Run File for the real band functions]
    ---------------------------------------------------------
%}


clear all;
close all;

%%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
D = 1; 
L1x = D;
L2x = 0;
L2y = D;
L1 = [L1x, 0];
L2 = [L2x,L2y];

%%% Radius of bubble
R=0.05;
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
c1 = D*[1/2,1/2];
c2 = D*[3/4,1/2];
N = 2; c = [c1;c2];

%%% Set the discretization or truncation parameters

N_multipole=3;
N_lattice=30;

N_XG = 60;
N_GM = floor(sqrt(2)*N_XG);
N_MX = N_XG;

%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;

%%% Plot the geometry
 figure, hold on
 t = linspace(0,2*pi);
 for i = -2:2
     for j = -2:2
         for n = 1:2
             plot(c(n,1)+R*cos(t)+i*L1(1)+j*L2(1), c(n,2)+R*sin(t)+i*L1(2)+j*L2(2),'k')  
             %text(c(n,1),c(n,2),num2str(n))
         end
     end
 end
 daspect([1 1 1])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set wave numbers ( k0: background, kb: bubble )
v0 = sqrt(kappa_0/rho_0);
vb = sqrt(kappa_b/rho_b);

k0 = 0.001;
JHdata = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);

alps = [(0:N_XG-1)*lengthXG/N_XG, (0:N_GM-1)*lengthGM/N_GM + lengthXG, (0:N_MX-1)*lengthMX/N_MX + lengthXG + lengthGM];
N_a = length(alps);
bandGM = zeros(N,N_GM);
bandMX = zeros(N,N_MX);
bandXG = zeros(N,N_XG);
    
for alp_ind = 1:N_XG
    %%% Set quasi-periodic parameter \alpha
    alp =  pointX + (alp_ind-1)/N_XG * (pointG-pointX);
    
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    ws = sort(vb*sqrt(delta*eig(C)./vol));
    bandXG(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_GM
    %%% Set quasi-periodic parameter \alpha
    alp =  pointG + (alp_ind-1)/N_GM * (pointM-pointG);
    
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    ws = sort(vb*sqrt(delta*eig(C)./vol));
    bandGM(:,alp_ind) = ws;
end
parfor alp_ind = 1:N_MX
    %%% Set quasi-periodic parameter \alpha
    alp =  pointM + (alp_ind-1)/N_MX * (pointX-pointM);
    
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    ws = sort(vb*sqrt(delta*eig(C)./vol));
    bandMX(:,alp_ind) = ws;
end
band0 = [bandXG, bandGM, bandMX];

% Plots the solution and saves good images
close all
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;

figure
hold on
w_real = sort(real(band0),1);

for I_w = 1:N
plot(alps,w_real(I_w,:),'b','LineWidth',lw)
end
x = [ 0*pi, lengthXG ,  lengthXG+lengthGM,      lengthXGM];
y= ['     X  ';'$\Gamma$';'    M   ';'   X    '];

axis([0,lengthXGM,0,max(max(w_real))])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 

set(gca,'fontsize', 15);
set(gca,'xgrid','on')

xlabel('Quasi-periodicity $\alpha$','interpreter','latex');
ylabel('Frequency $\omega$','interpreter','latex');
ax = gca;

set(gcf, 'Position', [0, 0, 600, 500])
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex')
l.FontSize = lfsz;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height*0.99];
