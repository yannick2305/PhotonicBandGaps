% computes the capacitance for a dimer


function matC = makeCR_2(k, z, R, alpha, beta, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice)   
                        
        N = 2; % two esonators inside of the unit cell

        matS = makeS(k,R,alpha,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice);

        % matR is now a 2x2 Block matrix with block entries that are
        % 2*N_SLP + 1 dimensional

        matR = makeR_2(k, z, R, alpha, beta, N_multi, N_lattice);

        matSR = (matS - matR)  ; % - matR)  ;    % Single Layer potential for beta neq 0
                                     % Minus sign because of the 2
                                     % conventions used (Habib vs Erik)

        % represent in Fourier Basis ---------------------------------
        n_range = -N_multi:N_multi;  % Define the range of the Fourier Basis
        coefficients_P = zeros(size(n_range));
        coefficients_M = zeros(size(n_range));

        % Compute the Fourier coefficients numerically
        %for i = 1:length(n_range)
        %    n = n_range(i);
        %    integrand = @(x) exp(R*beta(1)*cos(x) + R*beta(2)*sin(x)) .* exp(-1i*n*x);
        %    coefficients_P(i) = 1/(2*pi) * integral(integrand, 0, 2*pi);
        %end 

        %for i = 1:length(n_range)
        %    n = n_range(i);
        %    integrand = @(x) exp(-R*beta(1)*cos(x) - R*beta(2)*sin(x)) .* exp(-1i*n*x);
        %    coefficients_M(i) = 1/(2*pi) * integral(integrand, 0, 2*pi);
        %end 

         % Compute the Fourier coefficients using modified Bessel functions
        for i = 1:length(n_range)
            n = n_range(i);
            bI = besseli(n,R*norm(beta));
            coefficients_P(i) = exp(-1i*n*atan2(beta(2),beta(1)))*bI;
            coefficients_M(i) = exp(-1i*n*atan2(-beta(2),-beta(1)))*bI;
        end 
        %--------------------------------------------------------------
        NN = 2*N_multi + 1;
        MM = NN * N;

        matC = zeros(N,N);
        for j = 1:N
            phi_j = zeros(MM,1);
            phi_j((j-1)*NN+1 : j*NN) = coefficients_M; 
            psi_j = matSR\phi_j;
            for i = 1:N
                phi_i = zeros(MM,1);
                phi_i((i-1)*NN+1 : i*NN) = coefficients_P;
                matC(i,j) = -2*pi*R * phi_i' * psi_j;
            end
        end
end
