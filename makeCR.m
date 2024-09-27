
function matC = makeCR(k, R, alpha, beta, L1x, L2, d_zeta, JHdata, JHijdata, N, N_multi, N_lattice)   
        N = 1;
        matS = makeS(k,R,alpha,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice);
        matR = makeR(k, R, alpha, beta, N_multi, N_lattice);
        matSR = (-matR + matS)  ;    % Single Layer potential for beta neq 0
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

        phi_j = coefficients_P; 
        phi_j = transpose(phi_j);        
        psi_j = matSR\phi_j;                                % "Inverting" matSR
        matC = -2 * pi * R * dot(coefficients_M, psi_j);    % Compute capacitance matrix

end

