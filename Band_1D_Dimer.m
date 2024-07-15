% Spectral Plot for a Dimer

delta = 0.1;
s_1 = 0.8;
s_2 = 2;

% Define the Bandfunctions
f_1 = @(x) sqrt(delta) * sqrt((1/s_1 + 1/s_2) + sqrt(1/(s_1)^2 + 1/(s_2)^2 + 2/(s_1*s_2)*cos(x)));
f_2 = @(x) sqrt(delta) * sqrt((1/s_1 + 1/s_2) - sqrt(1/(s_1)^2 + 1/(s_2)^2 + 2/(s_1*s_2)*cos(x)));
f_3 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) - sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 - 2/(s_1*s_2)*cosh(x)))) );
f_4 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) + sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 - 2/(s_1*s_2)*cosh(x)))) );
f_5 = @(x) sqrt(delta) * sqrt(abs((1/s_1 + 1/s_2) + sqrt(abs(1/(s_1)^2 + 1/(s_2)^2 + 2/(s_1*s_2)*cosh(x)))) );

% Define the range for x
beta_bound = acosh((1/2) * (s_2^2 + s_1^2)/(s_1*s_2));
x_range = -pi:0.08:pi;
beta_range = -beta_bound:0.0098:beta_bound+0.001 ;
beta_range_2 = -4:0.08:4;

% Calculate the values of the function for each point in the grid
y_values_1 = real(f_1(x_range));
y_values_2 = real(f_2(x_range));
y_values_3 = real(f_3(beta_range));
y_values_4 = real(f_4(beta_range));
y_values_5 = real(f_5(beta_range_2));


% Create the initial plot for the functions
figure;
lw = 4.5; % Linewidth

plot(x_range, y_values_1, 'k', 'LineWidth', lw);
hold on;
pbaspect([8 6 1]); % Set the ratio to 10:6
plot(x_range, y_values_2, 'k', 'LineWidth', lw);
plot(beta_range, y_values_3, 'r', 'LineWidth', lw);
plot(beta_range, y_values_4, 'r', 'LineWidth', lw);
plot(beta_range_2, y_values_5, 'r', 'LineWidth', lw);

y_val = f_1(pi);
line([-4, 4], [y_val, y_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
b_val = f_2(pi);
line([-4, 4], [b_val, b_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
g_val = f_1(0);
line([-4, 4], [g_val, g_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);

xlabel('\alpha and \beta respectively', 'FontSize', 18);
ylabel('\omega^{\alpha, \beta}', 'FontSize', 18);
legend('Band function', '', 'Gap function', '', '', 'Location', 'northeast', 'FontSize', 20);
set(gca, 'FontSize', 22); % Adjust tick font size
saveas(gcf, 'Band_1D_Dimer.pdf');
hold off;
