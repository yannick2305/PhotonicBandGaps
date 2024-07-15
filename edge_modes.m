% edge modes


% Consider a chain of dimers with a defect in the middle
% We sudy the transfer matrix over the entire array 

clear all

nummber_range = 1000;

% Define range of k values

range = 1;
k_values =linspace(0, range, nummber_range);

% Initialize lists to store eigenvalues
eigenvalues_real = zeros(nummber_range, 2);
eigenvalues_imag = zeros(nummber_range, 2);
product_of_eigenvalues = zeros(nummber_range, 1);
transmission = zeros(nummber_range, 1);


l1 = 1;
l2 = 2;
def = 2;
n = 1.8;
del = 0.1;

format long;

% Compute eigenvalues for each k
for i = 1:nummber_range
    k = k_values(i);
    % Generate the matrix
    %A = Chain(k, l1, l2, def, del, n);
    %A = Trans(k, def, del, n);
    DL = Dimer(k, l1, l2, del, n);
    DR = DimerR(k, l1, l2, del, n);
    power = 2;

    A =  matrix_power(DL, power) * Trans(k, def, del, n) *  matrix_power(DR, power);
    % Compute eigenvalues
    eigenvalues = eig(A);

    transmission(i, :) = 1 / (1 + abs(A(1,2))^2);
    
    % Store real and imaginary parts separately
    eigenvalues_real(i, :) = real(eigenvalues);
    eigenvalues_imag(i, :) = imag(eigenvalues);
    
    % Compute product of eigenvalues and store
    product_of_eigenvalues(i) = prod(eigenvalues);

end


hold on;
ylim([-3, 3]);
pbaspect([8 6 1]); % Set the ratio to 10:6
scatter(k_values, product_of_eigenvalues, 5, 'black', 'filled', 'DisplayName', '', 'SizeData', 20);
scatter(0, 0, 5, 'red',  'filled', 'DisplayName', '', 'SizeData', 4);  % Just to get the legend right
scatter(0, 0, 5, 'blue', 'filled', 'DisplayName', '$', 'SizeData', 4);

scatter(repmat(k_values', 1, 2), eigenvalues_imag,      5, 'blue', 'filled', 'DisplayName', '$\Im(\lambda(k))$',   'SizeData', 10);
scatter(repmat(k_values', 1, 2), eigenvalues_real(:,1), 5, 'red',  'filled', 'DisplayName', '$\Re(\lambda_1(k))$', 'SizeData', 10);
scatter(repmat(k_values', 1, 2), eigenvalues_real(:,2), 5, 'red',  'filled', 'DisplayName', '$\Re(\lambda_2(k))$', 'SizeData', 10);

scatter(k_values, transmission,      5, 'green', 'filled', 'DisplayName', '$\Im(\lambda(k))$',   'SizeData', 20);

hold off;

xlabel('k');
ylabel('');
%legend('Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 16); % Adjust tick font size
box on;
saveas(gcf, 'Eigval_plot.pdf');

function T = TR(k, a, del, n)

    A = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];
    B = [1 + del/n, 1 - del/n; 1 - del/n, 1 + del/n];
    C = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];

    T =  0.5 * A * B * C;
end


function T = TL(k, a, del, n)
    
    A = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];
    B = [1 + n/del, 1 - n/del; 1 - n/del, 1 + n/del];
    C = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];

    T =  0.5 * A * B * C;
end


function T = Trans(k, l, del, n)
    T = TL(k, l, del, n) * TR(k, l, del, n);
end 

function D = Defect(k)
    D = TL(k) * TR(k);
end

function D = Dimer(k, l1, l2, del, n)
    D = Trans(k, l1, del, n) * Trans(k, l2, del, n);

end

function D = DimerR(k, l1, l2, del, n)
    D =  Trans(k, l2, del, n) * Trans(k, l1, del, n);

end


function T = Chain(k, l1, l2, def, del, n)
    
    T = Dimer(k, l1, l2, del, n) * Dimer(k, l1, l2, del, n) * Trans(k, def, del, n) * DimerR(k, l1, l2, del, n) * DimerR(k, l1, l2, del, n);

end


function result = matrix_power(M, n)
    % Check if M is a square matrix
    [rows, cols] = size(M);
    if rows ~= cols
        error('Input matrix must be square.');
    end
    
    % Initialize the result as the identity matrix
    result = eye(rows);
    
    % Compute M^n
    for i = 1:n
        result = result * M;
    end
end

