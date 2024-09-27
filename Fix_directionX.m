clear all
close all
% Suppose that beta has fixed, specified direction and alpha traverses BZ.
% Uses Muller's method to find zeros of imag(C).

%% GtX: branch of alpha from Gamma to X
z0s = [0]; % Precomputed initial guesses for R = 0.05 at X
N0 = length(z0s);
Na = 200; % Adjust as needed
alphas = linspace(pi*1e-5, pi-0.001, Na);

% Parameters for Muller's method
distTol = 5e-5;
fTol = 1e-8;
iterMax = 50;
e = 1e-4;

betas = zeros(Na,N0);
ws = zeros(Na,N0);
for I0 = 1:N0
    for Ia = 1:Na
        alp = alphas(Ia)*[0,1];
        func = @(bet) imag(my_function(alp,bet));
        if Ia == 1        
            z0 = z0s(I0);
        else
            z0 =  betas(Ia-1,I0);
        end
        betas(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        ws(Ia,I0) = my_function(alp,betas(Ia,I0));
    end
end

% Add branch that diverges in the middle
z0 = -4.24;
alphas15 = linspace(pi*1e-5, 0.344, Na);
betas15 = zeros(Na,1);
ws15 = zeros(Na,1);
for Ia = 1:Na
    alp = alphas15(Ia)*[0,1];
    func = @(bet) imag(my_function(alp,bet));
    if Ia > 1
        z0 =  betas15(Ia-1);
    end
    betas15(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    ws15(Ia) = my_function(alp,betas15(Ia));
end

% Add branch that diverges in the middle
z0 = -1.97739;
alphas2 = linspace(0.24, 0.344, Na);
betas2 = zeros(Na,1);
ws2 = zeros(Na,1);
betas2(1) = z0; 
for Ia = 1:Na
    alp = alphas2(Ia)*[0,1];
    func = @(bet) imag(my_function(alp,bet));
    if Ia > 1
        z0 =  betas2(Ia-1);
        betas2(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    end
    ws2(Ia) = my_function(alp,betas2(Ia));
end

%% Add branch for fixed alpha = 0
Nb = 200;
betas3 = linspace(-6,0,Nb);
ws3 = zeros(Nb,1);
alpha3 = 1e-5*[0,pi];
for Ib = 1:Nb
    ws3(Ib) = my_function(alpha3,betas3(Ib));
    if abs(imag(ws3(Ib))) > 5e-2
        ws3(Ib) = NaN;
    end
end
ws3(end-10)=ws3(end-11);
betas3(end-10) = 0;

%% Add branch for fixed alpha = [0,pi]
ws4 = zeros(Nb,1);
alpha4 = [0,pi];
betas4 = linspace(-9,0,Nb);
for Ib = 1:Nb
    ws4(Ib) = my_function(alpha4,betas4(Ib));
    if abs(imag(ws4(Ib))) > 1e-1
        ws4(Ib) = NaN;
    end
end

%% MtG: branch of alpha from M to X
z0s = [0]; % Precomputed initial guesses for R = 0.05 at M
N0 = length(z0s);
Na = 200; % Adjust as needed
alphas5 = linspace(0.005, pi-0.001, Na);

betas5 = zeros(Na,N0);
ws5 = zeros(Na,N0);
for I0 = 1:N0
    for Ia = 1:Na
        talp = alphas5(Ia);
        alp = [pi-talp,pi];
        func = @(bet) imag(my_function(alp,bet));
        if Ia == 1        
            z0 = z0s(I0);
        else
            z0 =  betas5(Ia-1,I0);
        end
        betas5(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        ws5(Ia,I0) = my_function(alp,betas5(Ia,I0));
    end
end

%% GtM: branch of alpha from Gamma to M
z0s = [0]; % Precomputed initial guesses for R = 0.05 at X
N0 = length(z0s);
Na = 200; % Adjust as needed
alphas6 = linspace(pi*1e-5, pi-0.001, Na);

betas6 = zeros(Na,N0);
ws6 = zeros(Na,N0);
for I0 = 1:N0
    for Ia = 1:Na
        alp = alphas6(Ia)*[1,1];
        func = @(bet) imag(my_function(alp,bet));
        if Ia == 1        
            z0 = z0s(I0);
        else
            z0 =  betas6(Ia-1,I0);
        end
        betas6(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        ws6(Ia,I0) = my_function(alp,betas6(Ia,I0));
    end
end

% Add branch that diverges in the middle
z0 = -4.13;
alphas7 = linspace(pi*1e-5, 0.234, Na);
betas7 = zeros(Na,1);
ws7 = zeros(Na,1);
for Ia = 1:Na
    alp = alphas7(Ia)*[1,1];
    func = @(bet) imag(my_function(alp,bet));
    if Ia > 1
        z0 =  betas7(Ia-1);
    end
    betas7(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    ws7(Ia) = my_function(alp,betas7(Ia));
end

% Add branch that diverges in the middle
z0 = -1.92765;
alphas8 = linspace(0.165, 0.234, Na);
betas8 = zeros(Na,1);
ws8 = zeros(Na,1);
betas8(1) = z0; 
for Ia = 1:Na
    alp = alphas8(Ia)*[1,1];
    func = @(bet) imag(my_function(alp,bet));
    if Ia > 1
        z0 =  betas8(Ia-1);
        betas8(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    end
    ws8(Ia) = my_function(alp,betas8(Ia));
end

%% Add branch for fixed alpha = [pi,pi]
ws9 = zeros(Nb,1);
alpha9 = [pi,pi];
betas9 = linspace(-10,0,Nb);
for Ib = 1:Nb
    ws9(Ib) = my_function(alpha9,betas9(Ib));
    if abs(imag(ws9(Ib))) > 1e-1
        ws9(Ib) = NaN;
    end
end
return

%%
%load plotXBZ
close all
figure
fsz = 16;
lfsz = 12;
alw = 0.75;
lw = 1.5;
msz = 3;

hold on

plot(3*pi-alphas15,ws15,'r','linewidth',lw) 
plot(betas15,ws15,'r','linewidth',lw) 

plot(3*pi-alphas2,ws2,'r','linewidth',lw) 
plot(betas2,ws2,'r','linewidth',lw) 

plot(betas3(1:end-10),real(ws3(1:end-10)),'r','linewidth',lw)
%plot(alpha3(2)*ones(size(ws3(1:end-10))),real(ws3(1:end-10)),'r','linewidth',lw)
plot(alpha3(1)*ones(size(ws3(1:end-10))),[1.05;ws3(2:end-10)],'r','linewidth',lw)
plot(betas4,real(ws4),'r','linewidth',lw)
plot(3*pi-alpha4(2)*ones(size(ws4)),real(ws4),'r','linewidth',lw)



plot(alphas7,ws7,'r','linewidth',lw) 
plot(betas7,ws7,'r','linewidth',lw) 

plot(alphas8,ws8,'r','linewidth',lw) 
plot(betas8,ws8,'r','linewidth',lw) 

plot(betas9,real(ws9),'r','linewidth',lw)
plot(alpha9(2)*ones(size(ws9)),real(ws9),'r','linewidth',lw)


plot(3*pi-alphas,ws(:,1),'k','linewidth',lw) 
plot(betas(:,1),ws(:,1),'k','linewidth',lw) 
plot(alphas5+pi,ws5(:,1),'k','linewidth',lw)
plot(betas5(:,1),ws5(:,1),'k','linewidth',lw) 
plot(alphas6,ws6(:,1),'k','linewidth',lw) 
plot(betas6(:,1),ws6(:,1),'k','linewidth',lw) 

ylim([0,1.05])
xlim([-6,3*pi])

ylabel('Frequency $\omega$','interpreter'   ,'latex');
xlabel('Complex part $\beta\qquad\qquad\quad\qquad$ Real part $\alpha\qquad\qquad\quad$','interpreter','latex');
xticks([ -6, -4, -2, 0, pi, 2*pi, 3*pi])
xticklabels({'$-6$','$-4$','$-2$','$\Gamma$','$M$','$X$','$\Gamma$'})
ax = gca;
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex')
%legend([b1(1), r1(1)],'Full solution','Perturbation theory','interpreter','latex','location','southeast','FontSize',lfsz)
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position', [0, 0, 700, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left*1 bottom ax_width*0.98 ax_height*0.98];

%print('XfixdirBZ','-depsc');

% Define the function f
function ws = my_function(alpha,tbet)
    beta = real(tbet);

    slopeb = 0;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % Set all the Paramters to the same values as in the RUN_band file
    k0 = 0.00001;
    N = 1;
    R = 0.05;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = [c1];
    N_lattice = 4;      % Use about 5
    N_multi = 4;          % Use about 5
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 1e-3;
    vb = 1;

    % Compute the Bandfunctions
    CR = makeCR(k0, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));

end