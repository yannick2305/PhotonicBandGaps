% clear all
% close all
% Suppose that alpha and beta have the same direction. Uses Muller's method
% to find zeros of imag(C).
% 
% Generate initial guesses for small alpha
num_points = 300; % Adjust as needed
betasearch  = linspace(-2*pi, 0, num_points);
%betasearch  = linspace(-1.91, -1.89, num_points);
talp = 0.5;
alp = [talp,0];

W = arrayfun(@(y) minfunci(alp,y,2), betasearch); 

figure
hold on
plot(betasearch,W,'.b')
plot(betasearch,0*W,'r')
%ylim([-0.001,0.001])
return

N = 2;
z0s = [-3.189]; % Precomputed initial guesses for R = 0.05 at M
%z0s = [0]; % Precomputed initial guesses for R = 0.05 at X
N0 = length(z0s);
Na = 200; % Adjust as needed
alphas = linspace(0.001, pi-0.001, Na);

% Parameters for Muller's method
distTol = 5e-5;
fTol = 1e-8;
iterMax = 50;
e = 1e-4;

betas = zeros(Na,N0);
ws = zeros(Na,N0);
for I0 = 1:N0
    for Ia = 1:Na
        alp = alphas(Ia)*[1,1];
        func = @(bet) minfunc(alp,bet);
        if Ia == 1        
            z0 = z0s(I0);
        else
            z0 =  betas(Ia-1,I0);
        end
        betas(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        w = wfunc(alp,betas(Ia,I0));
        [~,I] = min(abs(imag(w)));
        ws(Ia,I0) = w(I);
    end
end

%% Add branch that diverges in the middle
z0 = -1.80;
Na15 = 200;
alphas15 = linspace(0.325, 2.58258, Na15);
betas15 = zeros(Na15,1);
ws15 = zeros(Na15,1);
for Ia = 1:Na15
    alp = alphas15(Ia)*[1,1];
    func = @(bet) minfunci(alp,bet,2);
    if Ia > 1
        z0 =  betas15(Ia-1);
    end
    betas15(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    w = wfunc(alp,betas15(Ia));
    [~,I] = min(abs(imag(w)));
    ws15(Ia) = w(I);
end
betas15(end) = -1.5842;
ws15(end) = 0.717041;


%% Add branch that diverges in the middle
z0 = -1.75686;
alphas2 = linspace(1.85, pi-0.001, Na); 
betas2 = zeros(Na,1);
ws2 = zeros(Na,1);
for Ia = 1:Na
    alp = alphas2(Ia)*[1,1];
    func = @(bet) minfunci(alp,bet,1);
    if Ia > 1
        z0 =  betas2(Ia-1);
    end
    betas2(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    w = wfunc(alp,betas2(Ia));
    [~,I] = min(abs(imag(w)));
    ws2(Ia) = w(I);
end

z0 = -1.66656;
alphas25 = linspace(1.85, 2.58258, Na);  
betas25 = zeros(Na,1);
ws25 = zeros(Na,1);
for Ia = 1:Na
    alp = alphas25(Ia)*[1,1];
    func = @(bet) minfunci(alp,bet,2);
    if Ia > 1
        z0 =  betas25(Ia-1);
    end
    betas25(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    w = wfunc(alp,betas25(Ia));
    [~,I] = min(abs(imag(w)));
    ws25(Ia) = w(I);
end

%%
Na2p25 = round(Na/2);
z0s = [-1.75686, -1.66656];
alptemp = linspace(1.826,1.85,Na2p25); 
betas2p25 = zeros(Na2p25,2);
ws2p25 = zeros(Na2p25,2);
for I0 = 1:2
    for Ia = Na2p25:-1:1
        alp = alptemp(Ia)*[1,1];
        func = @(bet) minfunci(alp,bet,I0);
        if Ia < Na2p25
            z0 =  betas2p25(Ia+1,I0);
        end
        betas2p25(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        w = wfunc(alp,betas2p25(Ia,I0));
        [~,I] = min(abs(imag(w)));
        ws2p25(Ia,I0) = w(I);
    end
end
ws2p = [ws2p25(:,1);ws2];
ws25p = [ws2p25(:,2);ws25];
betas2p = [betas2p25(:,1); betas2];
betas25p = [betas2p25(:,2); betas25];
alphas2p = [alptemp, alphas2];
alphas25p = [alptemp, alphas25];
alphas2p = [1.826, alphas2p(7:end-2)];
ws2p = [0.68;ws2p(7:end-2)];
betas2p = [-1.73;betas2p(7:end-2)];
alphas25p = [1.826, alphas25p];
ws25p = [0.68;ws25p];
betas25p = [-1.73;betas25p];

%% Add branch for fixed alpha = 0
Nb = 200;
betas3 = linspace(-5.22087,0,Nb);
betas3 = linspace(-6,0,Nb);
ws3 = zeros(Nb,2);
alpha3 = 1e-10*[pi,pi];
for Ib = 1:Nb
    ws3(Ib,:) = wfunc(alpha3,betas3(Ib));
    if abs(imag(ws3(Ib,1))) > 5e-2
        ws3(Ib,1) = NaN;
    end
    if abs(imag(ws3(Ib,2))) > 5e-2
        ws3(Ib,2) = NaN;
    end
end
% ws3(end-10)=ws3(end-11);
% betas3(end-10) = 0;

%% Add branch for fixed alpha = [pi,pi]
ws4 = zeros(Nb,2);
alpha4 = [pi,pi];
betas4 = linspace(-10,0,Nb);
for Ib = 1:Nb
    ws4(Ib,:) = wfunc(alpha4,betas4(Ib));
    if abs(imag(ws4(Ib,1))) > 1e-1
        ws4(Ib,1) = NaN;
    end
    if abs(imag(ws4(Ib,2))) > 1e-1
        ws4(Ib,2) = NaN;
    end
end

%% Add branch for fixed beta = 0
ws5 = zeros(Nb,2);
alphas5 = linspace(pi*1e-5, pi-0.001, Na);
betas5 = 0*alphas5;
for Ia = 1:Na
    ws5(Ia,:) = wfunc(alphas5(Ia)*[1,1],0);
end

%% GtX: Path from Gamma to X
z0s = [-6.153, -4.6]; % -4.35
alphas6 = linspace(0.05, pi-0.001, Na); 
betas6 = zeros(Na,2);
ws6 = zeros(Na,2);
fTol = 1e-1;
for I0 = 1:2
    for Ia = 1:Na
        alp = alphas6(Ia)*[1,0];
        func = @(bet) minfunc(alp,bet);
        if Ia == 1        
            z0 = z0s(I0);
        else
            z0 =  betas6(Ia-1,I0);
        end
        betas6(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        Gtol = 0;
        if I0 == 2 && alphas6(Ia) < Gtol % Just cheating a bit by regularizing around alpha = Gamma         
            a = Gtol;
            f = @(x) -0.5*(2*x/a-1).*((2*x/a-1).^2 - 3);
            alp = alphas6(Ia)*( (1-f(alphas6(Ia)))*[1,1] + f(alphas6(Ia))*[1,0]);
        end
        w = wfunc(alp,betas6(Ia,I0));
        [~,I] = min(abs(imag(w)));
        ws6(Ia,I0) = w(I);
    end
end
fTol = 1e-8;

%%
z0 = -1.90264;
z0 = -2;
alphas65 = linspace(0.39, pi-0.001, Na);
alphas65 = linspace(0.5, pi-0.001, Na);
betas65 = zeros(Na,1);
ws65 = zeros(Na,1);
fTol = 1e-2;
distTol = 1e-2;
for Ia = 1:Na
    alp = alphas65(Ia)*[1,0];
    func = @(bet) minfunci(alp,bet,2);
    if Ia > 1
        z0 =  betas65(Ia-1);
    end
    betas65(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    w = wfunc(alp,betas65(Ia));
    [~,I] = min(abs(imag(w)));
    ws65(Ia) = w(I);
end
fTol = 1e-8;
distTol = 5e-5;

%% Add branch for fixed alpha = [0,pi]
Nb7 = 500;
betas7 = [linspace(-6,-3.4,Nb7/2), linspace(-1.25,0,Nb7/2)];
ws7 = zeros(Nb7,2);
alpha7 = [pi,0];
for Ib = 1:Nb7
    ws7(Ib,:) = wfunc(alpha7,betas7(Ib));
    if abs(imag(ws7(Ib,1))) > 1e-2
        ws7(Ib,1) = NaN;
    end
    if abs(imag(ws7(Ib,2))) > 1e-2
        ws7(Ib,2) = NaN;
    end
end

%% Add branch for fixed beta = 0
ws8 = zeros(Na,2);
alphas8 = linspace(pi*1e-5, pi-0.001, Na);
betas8 = 0*alphas8;
for Ia = 1:Na
    ws8(Ia,:) = wfunc(alphas8(Ia)*[1,0],0);
end

ws9 = zeros(Na,2);
alphas9 = linspace(pi*1e-5, pi-0.001, Na);
betas9 = 0*alphas9;
for Ia = 1:Na
    alp = [pi,pi-alphas9(Ia)];
    ws9(Ia,:) = wfunc(alp,0);
end

%% MtX branch from M to X
z0s = [-4.72815];%,-2.3956,-1.40794];
alphas10 = linspace(0.05, pi-0.001, Na); 
betas10 = zeros(Na,2);
ws10 = zeros(Na,2);
fTol = 1e-2;
distTol = 1e-2;
for I0 = 1:1
    for Ia = 1:Na
        alp = [pi, pi-alphas10(Ia)];
        func = @(bet) minfunc(alp,bet);
        if Ia == 1        
            z0 = z0s(I0);
        else
            z0 =  betas10(Ia-1,I0);
        end
        betas10(Ia,I0) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        w = wfunc(alp,betas10(Ia,I0));
        [~,I] = min(abs(imag(w)));
        ws10(Ia,I0) = w(I);
    end
end        
fTol = 1e-8;
distTol = 5e-5;

%% Touching up
betas65(end-2:end) = NaN;
ws65(end-2:end) = NaN;
ws6(end,2)=NaN;
betas6(end,2)=NaN;
ws6(end-2:end,1)=NaN;
betas6(end-2:end,1)=NaN;
betas10(end-2:end,1) = NaN;
betas10(end,2) = NaN;
betas10(end-1:end,3) = NaN;
ws10(end-2:end,1) = NaN;
ws10(end,2) = NaN;
ws10(end-1:end,3) = NaN;

%%
%close all
%load dimerBZ
figure
fsz = 16;
lfsz = 12;
alw = 0.75;
lw = 1.5;
msz = 3;

hold on
plot(alphas(1:end-1),ws(1:end-1),'r','linewidth',lw) 
plot(betas(1:end-1),ws(1:end-1),'r','linewidth',lw) 
plot(alphas15,ws15,'r','linewidth',lw) 
plot(betas15,ws15,'r','linewidth',lw) 
plot(alphas2p,ws2p,'r','linewidth',lw) 
plot(betas2p,ws2p,'r','linewidth',lw) 
plot(alphas25p,ws25p,'r','linewidth',lw) 
plot(betas25p,ws25p,'r','linewidth',lw) 

for I = 1:N
plot(betas3,real(ws3(:,I)),'r','linewidth',lw)
plot(alpha3(2)*ones(size(ws3(:,I))),real(ws3(:,I)),'r','linewidth',lw)
plot(betas4,real(ws4(:,I)),'r','linewidth',lw)
plot(alpha4(2)*ones(size(ws4(:,I))),real(ws4(:,I)),'r','linewidth',lw)
plot(betas5,real(ws5(:,I)),'k','linewidth',lw)
plot(alphas5,real(ws5(:,I)),'k','linewidth',lw)

plot(betas6(:,I),real(ws6(:,I)),'r','linewidth',lw)
plot(3*pi-alphas6,real(ws6(:,I)),'r','linewidth',lw)

plot(betas7,real(ws7(:,I)),'r','linewidth',lw)
plot(2*pi*ones(size(ws7(:,I))),real(ws7(:,I)),'r','linewidth',lw)
plot(betas8,real(ws8(:,I)),'k','linewidth',lw)
plot(3*pi-alphas8,real(ws8(:,I)),'k','linewidth',lw)
plot(betas9,real(ws9(:,I)),'k','linewidth',lw)
plot(pi+alphas9,real(ws9(:,I)),'k','linewidth',lw)
end

%plot(betas65,real(ws65),'r','linewidth',lw)
plot(3*pi-alphas65,real(ws65),'r','linewidth',lw)
for I = 1:0
plot(betas10(:,I),real(ws10(:,I)),'r','linewidth',lw)
plot(pi+alphas10,real(ws10(:,I)),'r','linewidth',lw)    
end

ylim([0,1.1])
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

%print('dimerBZ','-depsc');

%% 
figure
a = Gtol;
f = @(x) -0.5*(2*x/a-1).*((2*x/a-1).^2 - 3);
sx = linspace(0,a,100);
plot(sx,f(sx)) 


% Define the function f
function out = minfunc(alpha,tbet)
    beta = real(tbet);

    slopeb = 1;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % Set all the Paramters to the same values as in the RUN_band file
    %%% Set positions of two bubbles in unit cell
    k0 = 0.00001;
    N = 1;
    R = 0.05;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = D*[2/4,1/2];
    c2 = D*[3/4,1/2];
    N = 2; 
    c = [c1;c2];
    z = c2-c1;
    N_lattice = 8;      % Use about 5
    N_multi = 6;          % Use about 5
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
    CR = makeCR_2(k0, z, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));
    [~,I] = min(abs(imag(ws)));
    out = imag(ws(I));

end

function out = minfunci(alpha,tbet,i)
    beta = real(tbet);

    slopeb = 1;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % Set all the Paramters to the same values as in the RUN_band file
    %%% Set positions of two bubbles in unit cell
    k0 = 0.00001;
    N = 1;
    R = 0.05;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = D*[2/4,1/2];
    c2 = D*[3/4,1/2];
    N = 2; 
    c = [c1;c2];
    z = c2-c1;
    N_lattice = 8;      % Use about 5
    N_multi = 6;          % Use about 5
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
    CR = makeCR_2(k0, z, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));
    out = imag(ws(i));
end

function ws = wfunc(alpha,tbet)
    beta = real(tbet);

    slopeb = 1;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % Set all the Paramters to the same values as in the RUN_band file
    %%% Set positions of two bubbles in unit cell
    k0 = 0.00001;
    N = 1;
    R = 0.05;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = D*[2/4,1/2];
    c2 = D*[3/4,1/2];
    N = 2; 
    c = [c1;c2];
    z = c2-c1;
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
    CR = makeCR_2(k0, z, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));
end
