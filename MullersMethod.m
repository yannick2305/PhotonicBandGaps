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
