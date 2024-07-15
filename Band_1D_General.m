% Bandfunction in the General regime of a one-dimensional reosnator chain

n = 1.8;    % wave number coupling
a = 0.2;    % length of resonator
del = 1;    % contrast | Recover subwavelength regime (del=0.1) from the general regime
L = 1;      % length of the unit cell

num_points = 1003;   % number of scatter points

%-----------------------------------------------------------------------

k_values = linspace(0,5.5, num_points);

alpha_results = zeros(length(k_values), 2);
beta_results  = zeros(length(k_values), 2);
sub_res = zeros(length(k_values), 2);

for i = 1:length(k_values)
    alpha_results(i, :) = alpha(k_values(i), a, del, n, L);
    beta_results(i, :)  = beta(k_values(i), a, del, n, L);
end

figure;
pbaspect([8 6 1]); % Set the ratio to 10:6
scatter(k_values, real(alpha_results(:, 1)), 'red',  'filled', 'DisplayName', '\alpha(k)', 'SizeData', 30);
hold on;
scatter(k_values, real(beta_results(:, 1)) , 'black','filled', 'DisplayName', '\beta(k)', 'SizeData', 30);
%scatter(k_values, real(sub_res(:, 1)) , 'green','filled', 'DisplayName', '\beta(k)', 'SizeData', 30);
scatter(k_values, real(alpha_results(:, 2)), 'red',  'filled', 'SizeData', 30);
scatter(k_values, real(beta_results(:, 2)) , 'black','filled', 'SizeData', 30);
xlabel('k', 'FontSize', 18);
ylabel('' , 'FontSize', 18);
legend('\alpha(k)', '\beta(k)', 'Location', 'southwest', 'FontSize', 22);
set(gca, 'FontSize', 22); % Adjust tick font size
box on;
%saveas(gcf, 'General_Bandplot.pdf', 'pdf');


function T = T1(k, a, del, n) % Transfer over first boundary

    A = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];
    B = [1 + n/del, 1 - n/del; 1 - n/del, 1 + n/del];
    C = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];

    T =  0.5 * A * B * C;
end

function T = T2(k, a, del, n) % Transfer over second boundary

    A = [exp(-1i * n * k * a), 0; 0, exp(1i * n * k * a)];
    B = [1 + del/n, 1 - del/n; 1 - del/n, 1 + del/n];
    C = [exp(1i * k * a), 0; 0, exp(-1i * k * a)];

    T =  0.5 * A * B * C;
end

function F = FL(k, L)

    F = [exp(-1i * k * L), 0 ; 0, exp(1i * k * L)];
end 

function M = M(k, a, del, n, L)

    M = T1(k, a, del, n) * T2(k, a, del, n) * FL(k,L);
end

function b = beta(k, a, del, n, L) 

    eigenvalues = eig(M(k, a, del, n, L));
    b = [1 / L * log(abs(eigenvalues(1))), 1 / L * log(abs(eigenvalues(2)))];
end 

function a = alpha(k, a, del, n, L)

    eigenvalues = eig(M(k, a, del, n, L));
    a = [-1/L * angle(eigenvalues(1)), -1/ L * angle(eigenvalues(2))];
end
