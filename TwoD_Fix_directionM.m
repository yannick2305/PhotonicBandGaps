%{
    ---------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2024]
    Description:  [Plot the complex band structure (Experimental)]
    ---------------------------------------------------------
%}

% Code supporting the figures from arXiv:2409.14537
% -> More efficient implmentation can be found in TwoD_ComplexBands.m

% Suppose that beta has fixed, specified direction and alpha traverses BZ.
% Uses Muller's method to find zeros of imag(C).

clear all
close all

% --- Number of plotting points ---
    Na = 80; 
    Nb = 80;
    alphas = linspace(pi*1e-5, pi-0.001, Na);

% --- Parameters for Muller's method ---
    distTol = 5e-5;
    fTol    = 1e-8;
    iterMax = 50;
    e       = 1e-4;

%% --- Compute the spectral Bands ---

% --- Branch of alpha from Gamma to M ---
    z0s = [-5.19, 0]; 
    N0 = length(z0s);
    
    betas = zeros(Na,N0);
    ws = zeros(Na,N0);
    for I0 = 1:N0
        for Ia = 1:Na
            alp = alphas(Ia)*[1,1];
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

% --- Add branch that diverges in the middle ---
    z0 = -1.346;
    alphas2 = linspace(0.1699, pi-0.001, Na);
    betas2 = zeros(Na,1);
    ws2 = zeros(Na,1);
    betas2(1) = z0; 
    for Ia = 1:Na
        alp = alphas2(Ia)*[1,1];
        func = @(bet) imag(my_function(alp,bet));
        if Ia > 1
            z0 =  betas2(Ia-1);
            betas2(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        end
        ws2(Ia) = my_function(alp,betas2(Ia));
    end

% --- Add branch at alpha = 0 using Kummer's method ---
    betas3 = linspace(-6,0,Nb);
    ws3 = zeros(Nb,1);
    alpha3 = [0, 0];

    N_mul = 2;
    N_lat = 5;
    k0 = 0;
    R = 0.05; 
    vol = pi*R^2;
    delta = 1e-3;
    vb = 1;
    slope = 1;
   
    for i = 1:length(betas3)
        beta = betas3(i);            
        CR = makeCRKummer(k0, R, beta*[1,slope], N_mul, N_lat); 
        ws1 = abs(real(sort(vb * sqrt(abs(delta * eig(CR) ./ vol)))));      
    
        ws3(i) = ws1(1); 
    end

% --- Branch for fixed alpha = [pi,pi] ---
    ws4 = zeros(Nb,1);
    alpha4 = [pi,pi];
    betas4 = linspace(-3.2,0,Nb);
    for Ib = 1:Nb
        ws4(Ib) = my_function(alpha4,betas4(Ib));
        if abs(imag(ws4(Ib))) > 1e-1
            ws4(Ib) = NaN;
        end
    end

% --- Branch of alpha from M to X ---
    z0s = [-2.525, 0]; 
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

% --- Branch of alpha from X to Gamma ---
    z0s = [-5.08, 0]; 
    N0 = length(z0s);
    alphas6 = linspace(pi*1e-5, pi-0.001, Na);
    
    betas6 = zeros(Na,N0);
    ws6 = zeros(Na,N0);
    for I0 = 1:N0
        for Ia = 1:Na
            alp = [0,alphas6(Ia)];
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

% --- Add branch that diverges in the middle ---
    z0 = -1.39;
    alphas7 = linspace(0.3, pi-0.001, Na);
    betas7 = zeros(Na,1);
    ws7 = zeros(Na,1);
    for Ia = 1:Na
        alp = [0,alphas7(Ia)];
        func = @(bet) imag(my_function(alp,bet));
        if Ia > 1
            z0 =  betas7(Ia-1);
        end
        betas7(Ia) =  real(MullersMethod(func, z0-2*e, z0-e, z0, iterMax, distTol, fTol));
        ws7(Ia) = my_function(alp,betas7(Ia));
    end

% --- Add branch for transition to alpha = 0 ---
    betas8 = linspace(-6,0,Nb);
    ws8 = zeros(Nb,1);
    alpha8 = 1e-5*[0,pi];
    for Ib = 1:Nb
        ws8(Ib) = my_function(alpha8,betas8(Ib));
        if abs(imag(ws8(Ib))) > 5e-2
            ws8(Ib) = NaN;
        end
    end
    ws8(end-8)=ws8(end-9);
    betas8(end-8) = 0;

% --- Add branch for fixed alpha = [0,pi] ---
    ws9 = zeros(Nb,1);
    alpha9 = [0,pi];
    betas9 = linspace(-10,0,Nb);
    for Ib = 1:Nb
        ws9(Ib) = my_function(alpha9,betas9(Ib));
        if abs(imag(ws9(Ib))) > 1e-1
            ws9(Ib) = NaN;
        end
    end

%% --- Plot the spectral Bands ---
    close all
    figure
    fsz = 16;
    lfsz = 12;
    alw = 0.75;
    lw = 1.5;
    msz = 3;
    
    hold on
    plot(alphas,ws(:,1),'r','linewidth',lw) 
    plot(betas(:,1),ws(:,1),'r','linewidth',lw) 
    
    
    plot(alphas2,[ws2(1:end-3);ws2(end-3);ws2(end-3);ws2(end-3)],'r','linewidth',lw) 
    plot(betas2(1:end-1),ws2(1:end-1),'r','linewidth',lw) 
    
    plot(betas3,real(ws3),'r','linewidth',lw)
    plot(alpha3(1)*ones(size(ws3)),[1.05;ws3(2:end)],'r','linewidth',lw)
    plot(betas4,real(ws4),'r','linewidth',lw)
    plot(alpha4(1)*ones(size(ws4)),real(ws4),'r','linewidth',lw)
    
    plot([alphas5(1:end-14),alphas5(end)]+pi,[ws5(5,1);ws5(5,1);ws5(5,1);ws5(5,1);ws5(5:end-14,1);ws5(end-14,1)],'r','linewidth',lw) %check cont
    plot(betas5(2:end-14,1),ws5(2:end-14,1),'r','linewidth',lw)
    
    
    plot(3*pi-alphas6(1:end-3),ws6(1:end-3,1),'r','linewidth',lw)
    plot(betas6(1:end-3,1),ws6(1:end-3,1),'r','linewidth',lw)
    
    plot(3*pi-alphas7(1:end-3),ws7(1:end-3),'r','linewidth',lw)
    plot(betas7(1:end-3),ws7(1:end-3),'r','linewidth',lw)
    
    plot(betas9,real(ws9),'r','linewidth',lw)
    plot(alpha9(2)*ones(size(ws9))+pi,real(ws9),'r','linewidth',lw)
    
    
    plot(alphas,ws(:,2),'k','linewidth',lw) 
    plot(betas(:,2),ws(:,2),'k','linewidth',lw) 
    plot(3*pi-alphas6,ws6(:,2),'k','linewidth',lw)
    plot(betas6(:,2),ws6(:,2),'k','linewidth',lw)
    plot(alphas5+pi,ws5(:,2),'k','linewidth',lw)
    plot(betas5(:,2),ws5(:,2),'k','linewidth',lw) 
    
    
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
    
    %print('MfixdirBZ','-depsc');


%% --- Define the functions ---
% Define the function f
function ws = my_function(alpha,tbet)
    beta = real(tbet);

    slopeb = 1;

    alp = alpha;
    bet = [beta * slopeb, beta];

    % Set all the Paramters to the same values as in the RUN_band file
    k0 = 0.00001;
    N = 1;
    R = 0.05;     % R = 0.005 (Default setting)
    D = 1;         % D = 1 has to be true
    c1 = 1/2*D*[1,1];
    c = [c1];
    N_lattice = 4;       % Use about 5
    N_multi = 2;          % Use about 5
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
    CR = makeCRSlow(k0, R, alp, bet, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice);
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
        discriminant = sqrt(b^2 - 4 * a * c);
        
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
            %fprintf('Root found at x = %f after %d iterations.\n', real(root), i);
            return;
        end
        
        % Update points for the next iteration
        x0 = x1;
        x1 = x2;
        x2 = x3;
    end
    
    % If we reach the maximum iterations without convergence
    %warning('Maximum iterations reached without convergence.');
    root = x3;  % Return the best estimate of the root
end
