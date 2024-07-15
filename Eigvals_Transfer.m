% Define parameters
a = 1;
b = 2;
n = 1.8;
delta = 1; % has to be equal to 1 (General Regime)
l = 4;
n_0 = 1;

nummber_range = 600;

% Define range of k values
k_values =linspace(0, 2, nummber_range);

% Initialize lists to store eigenvalues
eigenvalues_real = zeros(nummber_range, 2);
eigenvalues_imag = zeros(nummber_range, 2);
product_of_eigenvalues = zeros(nummber_range, 1);

% Compute eigenvalues for each k
for i = 1:nummber_range
    k = k_values(i);
    % Generate the matrix
    %kp = k * n_0;
    A = matrix_element_2(k, n, a, b, delta);
    
    % Compute eigenvalues
    eigenvalues = eig(A);
    
    % Store real and imaginary parts separately
    eigenvalues_real(i, :) = real(eigenvalues);
    eigenvalues_imag(i, :) = imag(eigenvalues);
    
    % Compute product of eigenvalues and store
    product_of_eigenvalues(i) = prod(eigenvalues);

end

hold on;
pbaspect([8 6 1]); % Set the ratio to 10:6
scatter(k_values, product_of_eigenvalues, 5, 'black', 'filled', 'DisplayName', '$\lambda_1(k)\cdot\lambda_2(k)$', 'SizeData', 30);
scatter(0, 0, 5, 'red',  'filled', 'DisplayName', '$\Re(\lambda(k))$', 'SizeData', 4);  % Just to get the legend right
scatter(0, 0, 5, 'blue', 'filled', 'DisplayName', '$\Im(\lambda(k))$', 'SizeData', 4);

scatter(repmat(k_values', 1, 2), eigenvalues_real(:,1), 5, 'red',  'filled', 'DisplayName', '$\Re(\lambda_1(k))$', 'SizeData', 30);
scatter(repmat(k_values', 1, 2), eigenvalues_imag,      5, 'blue', 'filled', 'DisplayName', '$\Im(\lambda(k))$',   'SizeData', 30);
scatter(repmat(k_values', 1, 2), eigenvalues_real(:,2), 5, 'red',  'filled', 'DisplayName', '$\Re(\lambda_2(k))$', 'SizeData', 30);
hold off;

xlabel('k');
ylabel('');
legend('Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 16); % Adjust tick font size
box on;
saveas(gcf, 'Eigval_plot.pdf');

% Define matrix_element_2 function
function M = matrix_element_2(k, n, a, b, delta)
    im = 1i; % imaginary unit
    
    exp_term1 = exp(-2j * k * (b - a));
    %exp_term2 = exp(2j * b * k);
    %exp_term3 = exp(-2j * b * k);
    exp_term4 = exp(2j * k * (b - a));
    
    sin_term = sin(2 * a * k * n);
    cos_term = cos(2 * a * k * n);

    element1 = (-im * (delta^2 + n^2) * sin_term + 2 * delta * n * cos_term) / (2 * delta * n) * exp_term1;
    element2 = (im * (delta^2 - n^2)  * sin_term) / (2 * delta * n) * exp(2j * b * k);
    element3 = (-im * (delta^2 - n^2) * sin_term) / (2 * delta * n) * exp(-2j * b * k);
    element4 = ((delta^2 + n^2) * im  * sin_term + 2 * delta * n * cos_term) / (2 * delta * n) * exp_term4;
    
    M = [element1, element2; element3, element4];
end

