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
%   -> uncomment line 280
%
% Objective 3:
%   Specifically excite spectral bands using alpha-quasiperiodic sources
%   and overlay the result onto the spectral plot.
%   -> uncomment line 286

clear all;
close all;

% --- Set plotting parameters --- 
    Na = 50;  % Number of plotting points for real bandstructure.
    Nb = 50;  % Number of plotting points for complex bandstructure.

    % Parameters for Muller's method
    distTol = 5e-5;
    fTol    = 1e-8;
    iterMax = 50;
    e = 1e-4;


%% --- Compute the real bandstructure ---

% --- Branch of alpha from Gamma to X ---
    
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

% --- Branch of alpha from M to X ---
    
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

% --- Branch of alpha from Gamma to M ---
    
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


%% --- Compute the complex bandstructure ---

% --- Add branch at alpha = 0 using Kummer's method ---
    
    % Beta range
    betas3 = linspace(-6, 0, Nb); 
    ws3 = zeros(Nb,1);
    
    N_mul = 2;
    N_lat = 5;
    k0 = 0;
    R = 0.05; 
    vol = pi*R^2;
    delta = 1e-3;
    vb = 1;
    slope = 0;
   
    for i = 1:length(betas3)
        beta = betas3(i);            
        CR = makeCRKummer(k0, R, beta*[1,slope], N_mul, N_lat); 
        ws1 = abs(real(sort(vb * sqrt(abs(delta * eig(CR) ./ vol)))));      
    
        ws3(i) = ws1(1); 
    end
    
    
% --- Add branch for fixed alpha = [0,pi] ---
    ws4 = zeros(Nb,1);
    alpha4 = [0,pi];
    betas4 = linspace(-12, 0, Nb);
    for Ib = 1:Nb
        ws4(Ib) = (my_function(alpha4,betas4(Ib)));
        if abs(imag(ws4(Ib))) > 1e-2
            ws4(Ib) = NaN;
        end
    end


% --- Add branch for fixed alpha = [pi,pi] = M ---
    ws9 = zeros(Nb, 1);
    alpha9 = [pi, pi];
    betas9 = linspace(-12, 0, Nb);
    for Ib = 1:Nb
        ws9(Ib) = (my_function(alpha9, betas9(Ib)));
        if abs(imag(ws9(Ib))) > 1e-2
            ws9(Ib) = NaN;
        end
    end


%% --- Plotting the spectral Bands ---

% --- Precomputed values to overlay with spectral plot ---

    % --- Decay length of 2D defect eigenmode ---
    % (Precomputed values from: TwoD_DefectModes.m)
    decay_rate = [
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
       ];
    
    % --- alpha = [pi, pi] Quasiperiodic source ---
    % (Precomputed values from: TwoD_DefectMaster.m)
    quasi_source = [  
      -0.0004188   0.579973
      -0.9795922   0.588950
      -1.6763926   0.609973
      -2.2354314   0.639973
      -2.6216314   0.669973
      -2.9270593   0.699973
      -3.1861269   0.729973
      -3.4159571   0.759973
      -3.6266108   0.789973
      -3.8247867   0.819973
      -4.0154815   0.849973
      -4.2028857   0.879973
      -4.3909838   0.909973
      -4.5840896   0.939973
      -4.7874904   0.969973
      -5.0084023   0.999973
      -5.2574319   1.029973
      -5.5494919   1.059973
      ];

% --- Plotting the spectral Bands ---
    figure;
    fsz  = 19;
    lfsz = 14;
    alw  = 0.75;
    lw   = 1.5;
    msz  = 3;
    hold on

% --- Plot the complex bands ---
    
    % Band alpha = [0,0] (Kummer's method)
    plot(betas3, ws3,'r','linewidth',lw)
    line([0, 0]          , [0.578, 2.5], 'Color', 'r', 'LineStyle', '-', 'LineWidth', lw/2);
    line([9.4248, 9.4248], [0, 2.5],     'Color', 'r', 'LineStyle', '-', 'LineWidth', lw/2);
    
    % Band alpha = [0, pi]
    plot(betas4                        ,real(ws4),'r','linewidth',lw)
    plot(3*pi-alpha4(2)*ones(size(ws4)),real(ws4),'r','linewidth',lw/2)
    
    % Band alpha = [pi, pi]
    plot(betas9                   ,real(ws9),'r','linewidth',lw)
    plot(alpha9(2)*ones(size(ws9)),real(ws9),'r','linewidth',lw/2)

% --- Plot the real Band structure --- 

    % Path X -> G
    plot(3*pi-alphas,ws(:,1) ,'k','linewidth',lw)   
    plot(betas(:,1) ,ws(:,1) ,'k','linewidth',lw/2) 
    
    % Path M -> X
    plot(alphas5+pi ,ws5(:,1),'k','linewidth',lw)    
    plot(betas5(:,1),ws5(:,1),'k','linewidth',lw/2) 
    
    % Path G -> M
    plot(alphas6    ,ws6(:,1),'k','linewidth',lw)    
    plot(betas6(:,1),ws6(:,1),'k','linewidth',lw/2) 

% --- Uncomment below to get the desired result ---

    % --- Plot the decay as a function of the frequency ---
    constant = 0.11;
    %plot(decay_rate(:,1) + constant, decay_rate(:,2)    ,'x', 'Color', 'b', 'MarkerSize', 6,'LineWidth', lw );

    % Remark: The constant mereley serves as a normalisation and does not
    %         influence the exponential decay rate.

    % --- Excite specific bands using quasiperiodic sources ---
    %plot(quasi_source(:,1)*1.05 , quasi_source(:,2), 'x', 'Color', 'b', 'MarkerSize', 6,'LineWidth', lw);


    % --- Phase transition ---
    %{
    transp = 0.12;
    % --- Add Background colouring ---

    alM = 0.61;
    alG = 0.82;
    patch([-6 10 10 -6], [alG alG 2 2], [1, 0.647, 0]  , 'FaceAlpha', transp+0.2);
    patch([-6 10 10 -6], [alM alM alG alG],   'g', 'FaceAlpha', transp );
    patch([-6 10 10 -6], [0.579972 0.579972 alM alM], 'c', 'FaceAlpha', transp +0.35);

        % Plot each vertical line with different colors
    plot([0, 0]      , [0, 10], 'Color', [1, 0.647, 0], 'LineWidth', 4); 
    plot([pi, pi]    , [0, 10], 'Color', 'c',           'LineWidth', 4);
    plot([2*pi, 2*pi], [0, 10], 'Color', 'g',           'LineWidth', 4); 
    plot([3*pi, 3*pi], [0, 10], 'Color', [1, 0.647, 0], 'LineWidth', 4);

    fsf = 22 ;
    text(8, 0.593, '$\alpha = M$',      'Interpreter', 'latex', 'FontSize', fsf, 'Color', 'k');
    text(8, 0.71,  '$\alpha = X$',      'Interpreter', 'latex', 'FontSize', fsf, 'Color', 'k');
    text(8, 0.88,  '$\alpha = \Gamma$', 'Interpreter', 'latex', 'FontSize', fsf, 'Color', 'k');
    %}

% --- Figure and axis properties --- 
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
    set(gcf, 'Position', [0, 0, 700, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left*1 bottom ax_width*0.98 ax_height*0.98];

    % exportgraphics(gcf, 'XfixdirBZ.pdf'); % Save as a PDF file


%% --- Define the function and Muller's Method --- 

function ws = my_function(alpha,tbet)
    beta = real(tbet);

    slopeb = 0;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % --- Set the Parameters (same values as in RUN_band file) ---
    k0 = 0.001;
    N = 1;
    R = 0.05;           % Resonator size
    D = 1;        
    c1 = 1/2*D*[1,1];
    c = c1;
    N_lattice = 5;      % Use 5 (expectred truncation error of 1e-5)
    N_multi   = 2;      % Use 2 (expectred truncation error of 1e-7)
    d_zeta=makezetadata;
    JHdata = makeJHdata0(k0,R,N_multi);
    JHijdata = makeJHijexpdata(k0,c,N_multi);
    L1x = D;
    L2x = 0;
    L2y = D;
    L2 = [L2x,L2y];
    vol = pi*R^2;
    delta = 1e-3;
    vb = 1;

    % --- Compute the Bandfunctions --- 

        % Efficient tranformed lattice sum (Convergence O(N_lattice^{-3}))
        CR = makeCR(k0, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice); 
    
        % Non-efficient transformed method (Convergence O(N_lattice^{-3}))
        %CR = makeCRSlow(k0, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
    
        % Efficient direct computation     (Convergence O(N_lattice^{-1}))
        %CR = makeCRIntuitive(k0, R, alp, bet, N_multi, N_lattice);

        ws = sort(vb*sqrt(delta*eig(CR)./vol));

end
