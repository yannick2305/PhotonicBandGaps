%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Plot the complex band structure]
    ---------------------------------------------------------
%}

% Objective 1:
%   Generate a spectral plot of the complex bandstructure for beta along
%   a horizontal path, while alpha traverses the Brillouin zone. To compute 
%   the complex bandstructure, Muller's method is used to locate the zeros of 
%   the imaginary part of the capacitance matrix.
%
% Objective 2:
%   Overlay the spectral plot with the experimentally measured decay rates 
%   of the bandgap resonant modes. This overlay will demonstrate how the 
%   decay rates correspond to the complex bandstructure, illustrating the 
%   relationship between the observed decay and the complex band structure
%   of the system.
%
% Objective 3:
%   Specifically excite spectral bands using alpha-quasiperiodic sources
%   and overlay the result onto the spectral plot.

clear all;
close all;

% --- Set plotting parameters --- 
    Na = 50; % Number of plotting points
    Nb = 50;
    Nc = 50;
    % Parameters for Muller's method
    distTol = 5e-5;
    fTol    = 1e-8;
    iterMax = 50;
    e = 1e-4;

% --- Compute the spectral Bands ---

%% GtX: branch of alpha from Gamma to X

z0s = 0; % Precomputed initial guesses for R = 0.05 at X
N0 = length(z0s);

alphas = linspace(pi*1e-5, pi-0.001, Na);

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

%% Add branch that diverges in the middle

z0 = -4.24; % Precomputed initial guesses
alphas15 = linspace(pi*1e-4, 0.344, Na);
betas15  = zeros(Na,1);
ws15     = zeros(Na,1);

for Ia = 1:Na
    alp = alphas15(Ia)*[0,1];
    func = @(bet) imag(my_function(alp,bet));
    if Ia > 1
        z0 =  betas15(Ia-1);
    end
    betas15(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    ws15(Ia) = my_function(alp,betas15(Ia));
end

%% Add branch that diverges in the middle

z0 = -1.97739; % Precomputed initial guesses
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


%% Add branch at alpha = 0 using Kummer's method

% Beta range
%betas3 = linspace(-1.94, 0, Nc); 
%betas3 = linspace(-6, 0, Nc);
betas3 = linspace(-6, 0, Nc); 
ws3 = zeros(Nc,1);

N_mul = 2;
N_lat = 5;

 k0 = 0.00001;
    N = 1;
    R = 0.05; %R= 0.05;       % R = 0.005 (Default setting)
    D = 1;          % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = c1;
    N_lattice = 4;      % Use 4
    N_multi = 1;          % Use about 3- 4 and keep it fixed
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    %L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 1e-3;
    vb = 1;

    alpha = [0,0];

% Compute ws for each beta
parfor i = 1:length(betas3)
    beta = betas3(i);            
    CR = makeCRKummer(0, 0.05, beta*[1,0], N_mul, N_lat); 
    %CR = makeCR(0, R, alpha, beta*[1,0], L1x, L2, d_zeta, JHdata, JHijdata, N, N_mul, N_lat);
     ws = abs(real(sort(vb * sqrt(abs(delta * eig(CR) ./ vol)))));      

    ws3(i) = ws(1); % Store the smallest ws for plotting
end


%% Add branch for fixed alpha = 0, i.e. the unstable band
%{
betas3 = linspace(-6,0,Nc);
ws3 = zeros(Nc,1);
alpha3 = 1e-4*[0,1]; %1e-5*[0,pi];;
for Ib = 1:Nc
    ws3(Ib) = imag(my_function(alpha3,betas3(Ib)));
    if abs(imag(ws3(Ib))) > 3e-1
        ws3(Ib) = NaN;
    end
end
ws3(end-10)=ws3(end-11);
betas3(end-10) = 0;
%}

%% Add branch for fixed alpha = [0,pi]
ws4 = zeros(Nb,1);
alpha4 = [0,pi];
betas4 = linspace(-14,0,Nb);
for Ib = 1:Nb
    ws4(Ib) = (my_function(alpha4,betas4(Ib)));
    if abs(imag(ws4(Ib))) > 1e-1
        ws4(Ib) = NaN;
    end
end

%% MtG: branch of alpha from M to X

z0s = 0; % Precomputed initial guesses for R = 0.05 at M
N0 = length(z0s);
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

z0s = 0; % Precomputed initial guesses for R = 0.05 at X
N0 = length(z0s);
alphas6 = linspace(pi*1e-3, pi-0.001, Na);

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

%% Add transition branch that diverges in the middle
z0 = -4.13;
Nao = 200;
alphas7 = linspace(pi*1e-3, 0.2346, Nao);
betas7 = zeros(Nao,1);
ws7 = zeros(Nao,1);

% Parameters for Muller's method
distTol = 5e-6;
fTol = 1e-6;
iterMax = 150;
e = 1e-4;

for Ia = 1:Nao
    alp = alphas7(Ia)*[1,1];
    func = @(bet) imag(my_function(alp, bet));
    if Ia > 1
        z0 =  betas7(Ia-1);
    end
    betas7(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
    ws7(Ia) = my_function(alp,betas7(Ia));
end

%%
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


%% Add branch for fixed alpha = [pi,pi] = M
ws9 = zeros(Nb,1);
alpha9 = [pi,pi];
betas9 = linspace(-12,0,Nb); %linspace(-8*pi,0,Nb);
for Ib = 1:Nb
    ws9(Ib) = (my_function(alpha9,betas9(Ib)));
    if abs(imag(ws9(Ib))) > 1e-1
        ws9(Ib) = NaN;
    end
end



%% --- Plotting the spectral Bands ---

% --- Precomputed values to overlay with spectral plot---

    % Precomputed values from: DecayLength2D.m

    decay_rate_2 = [
      -0.43254   0.580
      -0.89855   0.585
      -1.12304   0.590
      -1.37771   0.600
      -1.67242   0.620
      -1.90621   0.640
      -2.12601   0.660
      -2.34337   0.680
      -2.56334   0.700
      -2.78602   0.720
      -2.99863   0.740
      -3.16206   0.760
      -3.27697   0.800
      -3.29682   0.840
      -3.27370   0.880
      -3.18269   0.920
      -3.05637   0.960
      -2.93697   1.000
      -2.83621   1.040
      -2.75303   1.080
      -2.68410   1.120
      -2.62640   1.160
      -2.57755   1.200
      -2.53575   1.240
      -2.49965   1.280
      -2.46819   1.320
      -2.44057   1.360
      -2.41615   1.400
      -2.39444   1.440
      -2.37502   1.480
      -2.35757   1.520
      -2.34182   1.560
      -2.32754   1.600
      -2.31455   1.640
      -2.30270   1.680
      -2.29184   1.720
      -2.28186   1.760
      -2.27267   1.800
      -2.26419   1.840
      -2.25634   1.880
      -2.24905   1.920
      -2.24228   1.960
      -2.23597   2.000
      -2.23008   2.040
      -2.22458   2.080
      -2.21943   2.120
      -2.21460   2.160
      -2.21006   2.200
      -2.20579   2.240
      -2.20177   2.280
      -2.19797   2.320
      -2.19439   2.360
      -2.19101   2.400
      -2.18780   2.440
      -2.18477   2.480];


    decay_new = [

   -0.25865   0.580
  -0.88923   0.620
  -1.49883   0.640
  -1.80710   0.660
  -2.06897   0.680
  -2.31781   0.700
  -2.56603   0.720
  -2.81920   0.740
  -3.07495   0.760
  -3.30833   0.780
  -3.39527   0.800
  -3.87634   0.840
  -3.65902   0.860
  -3.48027   0.880
  -3.32841   0.900
  -3.20052   0.920
  -3.09258   0.940
  -3.00074   0.960
  -2.92186   0.980
  -2.85344   1.000
  -2.79357   1.020
  -2.74076   1.040
  -2.69383   1.060
  -2.65187   1.080
  -2.61412   1.100
  -2.57999   1.120
  -2.54898   1.140
  -2.52070   1.160
  -2.49480   1.180
  -2.47099   1.200
  -2.44905   1.220
  -2.42876   1.240
  -2.40994   1.260
  -2.39245   1.280
  -2.37616   1.300
  -2.36094   1.320
  -2.34671   1.340
  -2.33336   1.360
  -2.32083   1.380
  -2.30904   1.400
  -2.29794   1.420
  -2.28746   1.440
  -2.27756   1.460
  -2.26819   1.480
  -2.25932   1.500
  -2.25090   1.520
  -2.24291   1.540
  -2.23531   1.560
  -2.22809   1.580
  -2.22120   1.600
  -2.21464   1.620
  -2.20837   1.640
  -2.20239   1.660
  -2.19667   1.680
  -2.19120   1.700
  -2.18596   1.720
  -2.18095   1.740
  -2.17613   1.760
  -2.17152   1.780
  -2.16709   1.800];
    
    % Precomputed values from: QuasiperiodicSource.m
    % alpha = [pi, pi] Quasiperiodic source
    quasi_source = [  
      -0.000418828534270   0.579973000000000
      -0.979592244356256   0.588950000000000
      -1.676392656278258   0.609973000000000
      -2.235431420971876   0.639973000000000
      -2.621631420683334   0.669973000000000
      -2.927059368783901   0.699973000000000
      -3.186126953933953   0.729973000000000
      -3.415957172537171   0.759973000000000
      -3.626610888988518   0.789973000000000
      -3.824786757593992   0.819973000000000
      -4.015481581392450   0.849973000000000
      -4.202885790892657   0.879973000000000
      -4.390983861837827   0.909973000000000
      -4.584089651303600   0.939973000000000
      -4.787490428840056   0.969973000000000
      -5.008402337016043   0.999973000000000
      -5.257431998378033   1.029973000000000
      -5.549491953381525   1.059973000000000
      ];
% alpha path

alpath = [
3.1416, 0.5800
3.1416, 0.5900
3.1416, 0.6000
3.1416, 0.6100
6.2832, 0.6200
6.2832, 0.6300
6.2832, 0.8000
6.2832, 0.8100
0.8251, 0.8200
0.7616, 0.8300
0.6664, 0.8400
0.5553, 0.8500
0.3649, 0.8600
0.0000, 0.8700
0.0000, 0.9100
0.0000, 1.1200
9.3930, 1.1300
9.2661, 1.1400
9.2344, 1.1500
9.1709, 1.1600
9.1392, 1.1700
9.1074, 1.1800
9.1074, 1.1900
9.0757, 1.2000
9.0440, 1.2100
9.0440, 1.2200
9.0122, 1.2300
9.0122, 1.2400
8.9805, 1.2500
8.9805, 1.2600
8.9805, 1.2700
8.9488, 1.2800
8.9488, 1.2900
8.9488, 1.3000
8.9170, 1.3100
8.9170, 1.3200
8.9170, 1.3300
8.9170, 1.3400
8.8853, 1.3500
8.8853, 1.3600
8.8853, 1.3700
8.8853, 1.3800
8.8853, 1.3900
8.8536, 1.4000
8.8536, 1.4100
8.8536, 1.4200
8.8536, 1.4300
8.8536, 1.4400
8.8536, 1.4500
8.8218, 1.4600
8.8218, 1.4700
8.8218, 1.4800
8.8218, 1.4900
8.8218, 1.5000
8.8218, 1.5100
8.8218, 1.5200
8.8218, 1.5300
8.8218, 1.5400
8.7901, 1.5500
8.7901, 1.6600
8.7901, 1.6700
8.7584, 1.6800
8.7584, 1.6900
8.7584, 1.7000
8.7584, 1.7100
8.7584, 1.7200
8.7584, 1.7300
8.7584, 1.7400
8.7584, 1.7500
8.7584, 1.7600
8.7584, 1.7700
8.7584, 1.7800
8.7584, 1.7900
8.7584, 1.8000
8.7584, 1.8100
8.7584, 1.8200
8.7584, 1.8300
8.7584, 1.8400
8.7584, 1.8500
8.7584, 1.8600
8.7266, 1.8700
8.7266, 2.2000];


quasi_test = [ -0.16956   0.580
  -1.02928   0.590
  -1.41381   0.600
  -1.69062   0.610
  -1.91176   0.620
  -2.09826   0.630
  -2.26092   0.640
  -2.40607   0.650
  -2.53783   0.660
  -2.65903   0.670
  -2.77177   0.680
  -2.87759   0.690
  -2.97764   0.700
  -3.07283   0.710
  -3.16381   0.720
  -3.25116   0.730
  -3.33534   0.740
  -3.41679   0.750
  -3.49594   0.760
  -3.57318   0.770
  -3.64889   0.780
  -3.72340   0.790
  -3.79695   0.800
  -3.86972   0.810
  -3.94178   0.820
  -4.01314   0.830
  -4.08384   0.840
  -4.15399   0.850
  -4.22389   0.860
  -4.29404   0.870
  -4.36498   0.880
  -4.43691   0.890
  -4.50922   0.900
  -4.58074   0.910
  -4.65390   0.920
  -4.75698   0.930];


% --- Plotting the spectral Bands ---
    figure;
    fsz  = 19;
    lfsz = 14;
    alw  = 0.75;
    lw   = 1.5;
    msz  = 3;
    hold on

% --- Plot the complex bands ---
    % Bands from 0 to [0,3] transition to singularity
    %plot(3*pi-alphas15,ws15,'r','linewidth',lw) 
    %plot(betas15,      ws15,'r','linewidth',lw) 
    %plot(betas2,        ws2,'r','linewidth',lw) 

    % Band alpha = [0,0]
    %plot(betas3(1:end-10),real(ws3(1:end-10)),'r','linewidth',lw)
    %plot(alpha3(2)*ones(size(ws3(1:end-10))),real(ws3(1:end-10)),'r','linewidth',lw)
    
    % Band alpha = [0,0] (Kummer's method)
    plot(betas3, ws3,'r','linewidth',lw)
    %plot(alpha3(1)*ones(size(ws3(1:end-10))),[1.05;ws3(2:end-10)],'r','linewidth',lw)
    line([0, 0]          , [0.578, 2.5], 'Color', 'r', 'LineStyle', '-', 'LineWidth', lw/1);
    line([9.4248, 9.4248], [0, 2.5],     'Color', 'r', 'LineStyle', '-', 'LineWidth', lw/1);
     

    %line([9.08078, 9.08078], [0.915,  2.5], 'Color', 'r', 'LineWidth', lw);
    
    % Band alpha = [0, pi]
    plot(betas4                        ,real(ws4),'r','linewidth',lw)
    plot(3*pi-alpha4(2)*ones(size(ws4)),real(ws4),'r','linewidth',lw/1)
    
    % Transition band to asymptote
    %plot(alphas7,ws7,'r','linewidth',lw/1) 
    %plot(betas7 ,ws7,'r','linewidth',lw) 
    
    % Band with singularity
    %plot(alphas8,ws8,'r','linewidth',lw/1) 
    %plot(betas8 ,ws8,'r','linewidth',lw) 
    
    % Band alpha = [pi, pi]
    plot(betas9                   ,real(ws9),'r','linewidth',lw)
    plot(alpha9(2)*ones(size(ws9)),real(ws9),'r','linewidth',lw/1)

% --- Plot the real Band structure --- 
    plot(3*pi-alphas,ws(:,1) ,'k','linewidth',lw)   % Path X -> G
    plot(betas(:,1) ,ws(:,1) ,'k','linewidth',lw/1) 
    plot(alphas5+pi ,ws5(:,1),'k','linewidth',lw)   % Path M -> X
    plot(betas5(:,1),ws5(:,1),'k','linewidth',lw/1) 
    plot(alphas6    ,ws6(:,1),'k','linewidth',lw)   % Path G ->M
    plot(betas6(:,1),ws6(:,1),'k','linewidth',lw/1) 

% --- Uncomment below to get the desired result ---
    % Plot the decay as a function of the frequency 
    %plot( decay_new(:,1)*1 , decay_new(:,2), 'x', 'Color', 'b', 'MarkerSize', 6,'LineWidth', 2*lw/3);
    %plot( alpath(:,1), alpath(:,2), '-', 'Color', 'b', 'MarkerSize', 6,'LineWidth', lw/2);
    plot( decay_rate_2(:,1) +0.13, decay_rate_2(:,2) ,'x', 'Color', 'b', 'MarkerSize', 6,'LineWidth', lw/2 );

    % Excite specific bands using quasiperiodic sources
    %plot( quasi_source(:,1)*1.04, quasi_source(:,2), 'x', 'Color', 'b', 'MarkerSize', 9,'LineWidth', lw/2);
    %plot( quasi_test(:,1), quasi_test(:,2), 'x', 'Color', 'b', 'MarkerSize', 6,'LineWidth', lw/2); % New values

% --- Figure and axis properties --- 
    % Adding labels, legend, ticks and figure properties
    ylim([0,1.8])     
    xlim([-6,3*pi])   
    
    ylabel('Frequency $\omega$','interpreter'   ,'latex');
    xlabel('Complex part $\beta\quad\qquad\qquad\qquad\qquad\qquad$ Real part $\alpha\qquad\qquad\quad$','interpreter','latex');
    xticks([ -6, -4, -2, 0, pi, 2*pi, 3*pi])
    xticklabels({'$-6$','$-4$','$-2$','$\Gamma$','$M$','$X$','$\Gamma$'})
    ax = gca;
    set(gca, 'FontSize', fsz, 'LineWidth', alw);
    set(gca,'TickLabelInterpreter','latex')
    %legend([b1(1), r1(1)],'Full solution','Perturbation theory','interpreter','latex','location','southeast','FontSize',lfsz)
    set(gca, 'FontSize', fsz, 'LineWidth', alw);
    set(gca,'TickLabelInterpreter','latex')
    %set(gcf, 'Position', [0, 0, 800, 500])
    set(gcf, 'Position', [0, 0, 700, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left*1 bottom ax_width*0.98 ax_height*0.98];
    
    %axis equal;

    exportgraphics(gcf, 'XfixdirBZ.pdf'); % Save as a PDF file

%print('XfixdirBZ','-depsc');

%% --- Define the function and Muller's Method --- 

function ws = my_function(alpha,tbet)
    beta = real(tbet);

    slopeb = 0;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % --- Set the Parameters (same values as in RUN_band file) ---
    k0 = 0.001;
    N = 1;
    R = 0.05; %R= 0.05;       % R = 0.005 (Default setting)
    D = 1;          % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = c1;
    N_lattice = 5;      % Use 4
    N_multi = 2;          % Use about 3 - 4 and keep it fixed
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    %L1 = [L1x, 0];
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 1e-3;
    vb = 1;

    % --- Compute the Bandfunctions --- 
    CR = makeCR(k0, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    ws = sort(vb*sqrt(delta*eig(CR)./vol));

end



function root = MullersMethod(func, z0_2e, z0_e, z0, iterMax, distTol, fTol)
    % MullersMethod implements Müller's method to find the root of a function.
    % 
    % Parameters:
    %   func     : Handle to the function whose root is being sought.
    %   z0_2e    : First initial guess (z0 - 2*e).
    %   z0_e     : Second initial guess (z0 - e).
    %   z0       : Third initial guess (z0).
    %   iterMax  : Maximum number of iterations.
    %   distTol  : Tolerance on the distance between consecutive iterates.
    %   fTol     : Tolerance on the function value at the root (|f(root)|).
    % 
    % Returns:
    %   root : The estimated root.
    
    % Initialize the initial guesses
    x0 = z0_2e;
    x1 = z0_e;
    x2 = z0;
    
    % Iterate using Müller's method
    for i = 1:iterMax
        % Calculate function values at the guesses
        f0 = func(x0);
        f1 = func(x1);
        f2 = func(x2);
        
        % Calculate h values
        h0 = x1 - x0;
        h1 = x2 - x1;
        
        % Calculate deltas
        delta0 = (f1 - f0) / h0;
        delta1 = (f2 - f1) / h1;
        
        % Calculate the second difference coefficient (a)
        d = (delta1 - delta0) / (h1 + h0);
        
        % Calculate coefficients
        a = d;
        b = delta1 + h1 * a;
        c = f2;
        
        % Calculate the discriminant
        discriminant = b^2 - 4 * a * c;
        
        % Handle the case when the discriminant is negative (complex roots)
        if discriminant < 0
            discriminant = sqrt(-discriminant) * 1i; % Take the imaginary part
        else
            discriminant = sqrt(discriminant); % Take the real part
        end
        
        % Choose the denominator for the quadratic formula (avoid cancellation)
        if abs(b + discriminant) > abs(b - discriminant)
            denominator = b + discriminant;
        else
            denominator = b - discriminant;
        end
        
        % Update the next guess (Müller's formula)
        dx = -2 * c / denominator;
        x3 = x2 + dx;
        
        % Check stopping criteria
        if abs(dx) < distTol || abs(func(x3)) < fTol
            root = x3;
            return;
        end
        
        % Update points for the next iteration
        x0 = x1;
        x1 = x2;
        x2 = x3;
    end
    
    % If we reach the maximum iterations without convergence
    warning('Maximum iterations reached without convergence.');
    root = x3;  % Return the best estimate of the root
end
